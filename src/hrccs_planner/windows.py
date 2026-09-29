"""Rank nights for dayside (emission) HRCCS.

For each target and night the usable time is: target above alt_min (and the
site's pointing limit, if any), Sun below the twilight altitude, planet on the
dayside half of the orbit (phase 0.25-0.75) and outside secondary eclipse.
For an eccentric orbit, the dayside and the phase curve use the orbital phase
(orbital angle from transit / 2 pi; ephemeris.eccentric_orbit), the eclipse is
excluded around its own time (ephemeris.eclipse_phase), and the planet velocity
is the Keplerian curve.

The planet signal in each exposure is scaled by a Lambertian phase curve f
(ephemeris.lambert_phase), and each exposure can be weighted by an airmass term
w(X) = 10^(-0.4 k (X-1)) X^(-p) (extinction, and slit losses for seeing wider
than the slit). A matched-filter sum over exposures gives

    S/N ~ sqrt(sum(f^2 w dt) * delta_v) * W

with delta_v the range of planet velocity covered and W an optional per-target
weight: 10^(-0.2 (Kmag - 6)) for photon-noise-limited stars, times the planet's
dayside contrast. S/N is normalized to the best window (overall, or per target).
With max_hours, each window is limited to the contiguous stretch of at most that
length that maximizes sum(f^2 w dt) * delta_v. See docs/method.md.
"""
import csv
import os
import re
import sys
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .ephemeris import eccentric_orbit, eclipse_phase, lambert_phase, orbital_phase, planet_rv_circular
from .visibility import airmass as airmass_of
from .visibility import moon_too_close, night_sky


@dataclass
class WindowOptions:
    duration: float = 0.            # eclipse duration [h], excluded around the eclipse (phase 0.5 if circular)
    alt_min: float = 30.            # minimum target altitude [deg] (airmass 2)
    twilight: float = -18.          # Sun altitude defining night [deg]
    delta_v_min: float = 30.        # minimum planet velocity change [km/s]
    bright_threshold: float = 0.5   # Moon illumination above which a night is bright time
    moon_sep_min: float = 10.       # minimum target-Moon separation [deg]
    inclination: float = 90.        # [deg]; Kp = Kp_max sin(i)
    weight_kmag: bool = False
    contrast_column: Optional[str] = None   # target-list column with the contrast weight
    max_hours: Optional[float] = None
    airmass_k: float = 0.           # extinction [mag/airmass]
    seeing_exponent: float = 0.
    normalize_per_target: bool = False
    snr_min: float = 0.3
    sort: str = 'moon'              # 'moon': bright time first, then S/N; 'snr': S/N only
    top: Optional[int] = None
    samples: int = 1441             # time samples per night (1441: 1 minute)


def best_subwindow(usable, w2, vpl, nmax):
    """Mask of the stretch of at most nmax consecutive samples that maximizes
    sum(w2 over usable samples) * (planet velocity range over usable samples),
    where w2 is the per-sample S/N^2 weight (phase curve f^2 x airmass weight)."""
    idx = np.flatnonzero(usable)
    i0, i1 = idx[0], idx[-1] + 1
    if i1 - i0 <= nmax:
        return usable
    u = usable[i0:i1]
    f2 = np.where(u, w2[i0:i1], 0.)
    v = np.where(u, vpl[i0:i1], np.nan)
    eff = sliding_window_view(f2, nmax).sum(axis=1)
    vw = sliding_window_view(v, nmax)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        dv = np.nan_to_num(np.nanmax(vw, axis=1) - np.nanmin(vw, axis=1))
    k = np.argmax(eff * dv)
    out = np.zeros_like(usable)
    out[i0 + k:i0 + k + nmax] = u[k:k + nmax]
    return out


def target_weight(target, opts):
    """Per-target S/N weight (NaN if a required value is missing)."""
    w = 1.
    if opts.weight_kmag:
        if not np.isfinite(target.kmag):
            return np.nan
        w *= 10 ** (-0.2 * (target.kmag - 6.))
    if opts.contrast_column:
        try:
            w *= float(target.extra[opts.contrast_column])
        except (KeyError, TypeError, ValueError):
            return np.nan
    return w


def rank_windows(targets, site, dates, opts=None, log=print):
    """Best window per target and night.

    Returns (rows, n_all): the windows passing opts.snr_min, sorted, as dicts
    (keys: date, name, bright, illum, hrs, eff_hrs, fpl, ut0, ut1, X_mean,
    X_max, ph0, ph1, v0, v1, dv, crosses_eclipse, weight, snr), and the number
    of windows before the S/N cut.
    """
    opts = opts or WindowOptions()
    if opts.contrast_column is not None:
        have = any(opts.contrast_column in t.extra for t in targets)
        if not have:
            raise ValueError(f"no '{opts.contrast_column}' column in the target list")
    usable_targets = []
    skipped = []
    for t in targets:
        if not np.isfinite(target_weight(t, opts)):
            skipped.append(t.name)
            continue
        usable_targets.append(t)
    if skipped:
        log('Skipping (missing Kmag or contrast): ' + ', '.join(skipped))
    coords = {t.name: t.coord for t in usable_targets}
    sini = np.sin(np.radians(opts.inclination))
    grid = site.night_hours(opts.samples)
    dt_hr = grid[1] - grid[0]

    rows = []
    for date in dates:
        night = night_sky(site, date, opts.samples)
        hours = night.hours
        dark = night.dark(opts.twilight)
        illum = night.moon_illumination(opts.twilight)
        bright = illum >= opts.bright_threshold
        for t in usable_targets:
            weight = target_weight(t, opts)
            kp = t.kp_max * sini
            ecl_halfwidth = 0.5 * opts.duration / 24. / t.period
            altaz = coords[t.name].transform_to(night.frame)
            alt = altaz.alt.deg
            X = airmass_of(alt)
            wam = 10 ** (-0.4 * opts.airmass_k * (X - 1.)) * X ** (-opts.seeing_exponent)
            if moon_too_close(altaz, night, opts.moon_sep_min):
                continue
            phase = orbital_phase(night.times.jd, t.t0, t.period)
            if t.eccentric:
                ### Orbital phase for the geometry; the eclipse at its own time
                uphase, _, rv = eccentric_orbit(phase, t.e, t.omega)
                vpl = kp * rv
                phase = uphase % 1
                dphase = (orbital_phase(night.times.jd, t.t0, t.period) - eclipse_phase(t.e, t.omega) + 0.5) % 1 - 0.5
                in_eclipse = np.abs(dphase) < ecl_halfwidth
            else:
                vpl = planet_rv_circular(phase, kp)
                in_eclipse = np.abs(phase - 0.5) < ecl_halfwidth
            fpl = lambert_phase(phase, opts.inclination)
            above = alt > opts.alt_min
            if site.min_altitude is not None:
                above &= alt > site.min_altitude(altaz.az.deg)
            dayside = (phase > 0.25) & (phase < 0.75)
            usable = dark & above & dayside & ~in_eclipse
            if not np.any(usable):
                continue
            if opts.max_hours is not None:
                usable = best_subwindow(usable, fpl ** 2 * wam, vpl, int(round(opts.max_hours / dt_hr)) + 1)
            eff_hrs = np.sum(fpl[usable] ** 2 * wam[usable]) * dt_hr
            delta_v = np.max(vpl[usable]) - np.min(vpl[usable])
            if delta_v < opts.delta_v_min:
                continue
            ut = hours[usable]
            span = (hours >= ut[0]) & (hours <= ut[-1])
            rows.append({'date': date, 'name': t.name, 'bright': bright, 'illum': illum,
                         'hrs': np.sum(usable) * dt_hr, 'eff_hrs': eff_hrs, 'fpl': np.mean(fpl[usable]),
                         'ut0': ut[0], 'ut1': ut[-1],
                         'X_mean': np.mean(X[usable]), 'X_max': np.max(X[usable]),
                         'ph0': phase[usable][0], 'ph1': phase[usable][-1],
                         'v0': vpl[usable][0], 'v1': vpl[usable][-1], 'dv': delta_v,
                         'crosses_eclipse': np.any(in_eclipse & dark & above & span), 'weight': weight})

    for r in rows:
        r['snr'] = np.sqrt(r['eff_hrs'] * r['dv']) * r['weight']
    if opts.normalize_per_target:
        best = {}
        for r in rows:
            best[r['name']] = max(best.get(r['name'], 0.), r['snr'])
        for r in rows:
            r['snr'] /= best[r['name']]
    else:
        best = max([r['snr'] for r in rows]) if rows else 1.
        for r in rows:
            r['snr'] /= best
    n_all = len(rows)
    rows = [r for r in rows if r['snr'] >= opts.snr_min]
    if opts.sort == 'snr':
        rows.sort(key=lambda r: -r['snr'])
    else:
        rows.sort(key=lambda r: (not r['bright'], -r['snr']))
    if opts.top is not None:
        rows = rows[:opts.top]
    return rows, n_all


CSV_HEADER = ['rank', 'hst_date', 'target', 'moon', 'illum', 'usable_hrs', 'ut_start', 'ut_end',
              'phase_start', 'phase_end', 'v_start', 'v_end', 'delta_v', 'f_pl', 'rel_snr',
              'crosses_eclipse', 'airmass_mean', 'airmass_max']
# 'hst_date' is the local date on which the night starts, at any site (the name
# is kept for compatibility with existing readers, e.g. KPIC-Synthgen).


def table_lines(rows, first):
    """Printable table of ranked windows, preceded by the line `first`."""
    lines = [first,
             f"{'rank':>4} {'Date':10} {'target':14} {'moon':6} {'illum':>5} {'hrs':>5} {'UT start':>8} "
             f"{'UT end':>7} {'ph0':>5} {'ph1':>5} {'v0':>6} {'v1':>6} {'dv':>5} {'f_pl':>5} {'X':>4} "
             f"{'Xmax':>4} {'S/N':>5} {'ecl':>4}"]
    for i, r in enumerate(rows):
        lines.append(
            f"{i + 1:4d} {r['date']:10} {r['name'][:14]:14} {'bright' if r['bright'] else 'dark':6} "
            f"{r['illum']:5.2f} {r['hrs']:5.2f} {r['ut0']:8.2f} {r['ut1']:7.2f} {r['ph0']:5.2f} "
            f"{r['ph1']:5.2f} {r['v0']:6.0f} {r['v1']:6.0f} {r['dv']:5.0f} {r['fpl']:5.2f} "
            f"{r['X_mean']:4.2f} {r['X_max']:4.2f} {r['snr']:5.2f} {'Y' if r['crosses_eclipse'] else '':>4}")
    return lines


def write_windows(path, rows, first, command=None):
    """Ranked windows to `path` (CSV) and the table, with the command line, to .txt."""
    command = command if command is not None else ' '.join(sys.argv)
    with open(os.path.splitext(path)[0] + '.txt', 'w') as fh:
        fh.write(command + '\n' + '\n'.join(table_lines(rows, first)) + '\n')
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(CSV_HEADER)
        for i, r in enumerate(rows):
            w.writerow([i + 1, r['date'], r['name'], 'bright' if r['bright'] else 'dark', f"{r['illum']:.2f}",
                        f"{r['hrs']:.2f}", f"{r['ut0']:.2f}", f"{r['ut1']:.2f}", f"{r['ph0']:.3f}",
                        f"{r['ph1']:.3f}", f"{r['v0']:.1f}", f"{r['v1']:.1f}", f"{r['dv']:.1f}",
                        f"{r['fpl']:.3f}", f"{r['snr']:.3f}", 'Y' if r['crosses_eclipse'] else 'N',
                        f"{r['X_mean']:.3f}", f"{r['X_max']:.3f}"])
    return [path, os.path.splitext(path)[0] + '.txt']


def gemini_timing_windows(rows, round_minutes=15):
    """Gemini PIT Scheduling-field text: one block per target,

        TW for <target>
        --------------------------------------------
        YYYY-MM-DD HH:MM:SS H:MM        (UT start, duration)

    with windows in time order. With round_minutes > 0 the start is rounded down
    and the end up to a multiple of round_minutes, so each window covers the
    whole usable time."""
    from datetime import datetime, timedelta
    blocks = []
    for name in dict.fromkeys(r['name'] for r in rows):
        lines = []
        for r in sorted((r for r in rows if r['name'] == name), key=lambda r: (r['date'], r['ut0'])):
            base = datetime.fromisoformat(r['date']) + timedelta(hours=23, minutes=59, seconds=59)
            t0 = base + timedelta(hours=float(r['ut0']))
            t1 = base + timedelta(hours=float(r['ut1']))
            # nearest minute (the time grid is offset by 1 s from whole minutes)
            t0, t1 = [(t + timedelta(seconds=30)).replace(second=0, microsecond=0) for t in (t0, t1)]
            if round_minutes:
                step = timedelta(minutes=round_minutes)
                t0 -= (t0 - t0.replace(hour=0, minute=0)) % step
                rem = (t1 - t1.replace(hour=0, minute=0)) % step
                t1 += (step - rem) if rem else timedelta(0)
            minutes = int(round((t1 - t0).total_seconds() / 60.))
            lines.append(f'{t0:%Y-%m-%d %H:%M:%S} {minutes // 60}:{minutes % 60:02d}')
        blocks.append('\n'.join([f'TW for {name}', '-' * 44] + lines))
    return '\n\n'.join(blocks) + '\n'


def split_path(path, name):
    stem, ext = os.path.splitext(path)
    return f"{stem}_{re.sub(r'[^A-Za-z0-9]', '', name)}{ext}"


def read_windows(path):
    """Rows of a ranked-window CSV as dicts of strings."""
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh))
