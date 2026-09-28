"""Rank transit windows for transmission HRCCS.

A transit window is the transit itself (first to fourth contact, T14 centred on
T_c = T0 + n P) plus an out-of-transit baseline of `baseline` hours in total:

    baseline_mode = 'split': half the baseline before the transit and half after;
    baseline_mode = 'any':   any division before/after, including all of it on one
                             side. The most even division that fits is used, and
                             among equally even ones the lowest mean airmass.

The whole window (transit and baseline) must be at night (Sun below the twilight
altitude) with the target above alt_min (airmass 2 for 30 deg) and the site's
pointing limit, if any. The Moon must stay at least moon_sep_min from the target.

The relative S/N counts only the in-transit time, weighted by airmass:

    S/N ~ sqrt(sum over in-transit exposures of w(X) dt) * W

with w(X) = 10^(-0.4 k (X-1)) X^(-p) as for emission windows, and the optional
per-target weight W (Kmag, a contrast/depth column). No phase curve or planet
velocity term is used. Transit times use T0 and the period only, so eccentric
orbits are fine as long as T0 is a transit time.
"""
import csv
import os
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .visibility import airmass as airmass_of
from .visibility import moon_too_close, night_sky
from .windows import target_weight


@dataclass
class TransitOptions:
    duration: Optional[float] = None       # T14 [h] for all targets
    duration_column: Optional[str] = None  # target-list column with T14 [h] (overrides `duration`)
    baseline: float = 2.                   # total out-of-transit baseline [h]
    baseline_mode: str = 'split'           # 'split' (half before, half after) or 'any'
    alt_min: float = 30.
    twilight: float = -18.
    bright_threshold: float = 0.5
    moon_sep_min: float = 10.
    weight_kmag: bool = False
    contrast_column: Optional[str] = None
    airmass_k: float = 0.
    seeing_exponent: float = 0.
    normalize_per_target: bool = False
    snr_min: float = 0.3
    sort: str = 'moon'
    top: Optional[int] = None
    samples: int = 1441


def _all_ok(ok, hours, a, b):
    """True if [a, b] (UT hours) lies inside the grid and every sample in it is ok."""
    if a < hours[0] or b > hours[-1]:
        return False
    inside = (hours >= a) & (hours <= b)
    return bool(np.any(inside)) and bool(np.all(ok[inside]))


def choose_baseline(ok, hours, t1, t4, baseline, mode, X=None):
    """(pre, post) baseline hours around the transit [t1, t4] such that the whole
    window [t1 - pre, t4 + post] is usable, or None if impossible.

    mode 'split': pre = post = baseline / 2.
    mode 'any': pre + post = baseline, most even division first; ties go to the
    lowest mean airmass X over the window (if given)."""
    if mode == 'split':
        pre = post = baseline / 2.
        return (pre, post) if _all_ok(ok, hours, t1 - pre, t4 + post) else None
    if mode != 'any':
        raise ValueError(f"baseline_mode must be 'split' or 'any', not '{mode}'")
    step = hours[1] - hours[0]
    pres = np.unique(np.clip(np.r_[np.arange(0., baseline, step), baseline], 0., baseline))
    best, key = None, None
    for pre in pres:
        post = baseline - pre
        a, b = t1 - pre, t4 + post
        if not _all_ok(ok, hours, a, b):
            continue
        xm = float(np.mean(X[(hours >= a) & (hours <= b)])) if X is not None else 0.
        k = (-min(pre, post), xm)
        if key is None or k < key:
            best, key = (float(pre), float(post)), k
    return best


def transit_duration(target, opts):
    if opts.duration_column:
        try:
            return float(target.extra[opts.duration_column])
        except (KeyError, TypeError, ValueError):
            pass
    return opts.duration


def rank_transits(targets, site, dates, opts=None, log=print):
    """Transit windows per target and night, ranked by relative S/N.

    Returns (rows, n_all). Row keys: date, name, bright, illum, tc (UT h),
    t14, pre, post, ut0, ut1 (window, UT h), hrs, X_mean and X_in_max (in
    transit), X_max (whole window), weight, snr."""
    opts = opts or TransitOptions()
    use = []
    for t in targets:
        if t.transiting is False:
            log(f'Skipping {t.name} - not transiting (Transiting? = N)')
            continue
        d = transit_duration(t, opts)
        if d is None or not np.isfinite(d) or d <= 0:
            log(f'Skipping {t.name} - no transit duration (give --duration or a --duration-column value)')
            continue
        if not np.isfinite(target_weight(t, opts)):
            log(f'Skipping {t.name} - missing Kmag or contrast')
            continue
        use.append((t, d))
    coords = {t.name: t.coord for t, _ in use}
    grid = site.night_hours(opts.samples)
    dt_hr = grid[1] - grid[0]

    rows = []
    for date in dates:
        night = night_sky(site, date, opts.samples)
        hours, jd = night.hours, night.times.jd
        dark = night.dark(opts.twilight)
        illum = night.moon_illumination(opts.twilight)
        bright = illum >= opts.bright_threshold
        for t, d in use:
            altaz = coords[t.name].transform_to(night.frame)
            if moon_too_close(altaz, night, opts.moon_sep_min):
                continue
            alt = altaz.alt.deg
            X = airmass_of(alt)
            wam = 10 ** (-0.4 * opts.airmass_k * (X - 1.)) * X ** (-opts.seeing_exponent)
            above = alt > opts.alt_min
            if site.min_altitude is not None:
                above &= alt > site.min_altitude(altaz.az.deg)
            ok = dark & above
            if not np.any(ok):
                continue
            for n in range(int(np.floor((jd[0] - t.t0) / t.period)), int(np.ceil((jd[-1] - t.t0) / t.period)) + 1):
                tc = hours[0] + (t.t0 + n * t.period - jd[0]) * 24.
                t1, t4 = tc - d / 2., tc + d / 2.
                if t1 < hours[0] or t4 > hours[-1]:
                    continue
                choice = choose_baseline(ok, hours, t1, t4, opts.baseline, opts.baseline_mode, X)
                if choice is None:
                    continue
                pre, post = choice
                intr = (hours >= t1) & (hours <= t4)
                win = (hours >= t1 - pre) & (hours <= t4 + post)
                rows.append({'date': date, 'name': t.name, 'bright': bright, 'illum': illum, 'tc': tc, 't14': d,
                             'pre': pre, 'post': post, 'ut0': t1 - pre, 'ut1': t4 + post, 'hrs': d + pre + post,
                             'eff_hrs': np.sum(wam[intr]) * dt_hr, 'X_mean': np.mean(X[intr]),
                             'X_in_max': np.max(X[intr]), 'X_max': np.max(X[win]),
                             'weight': target_weight(t, opts)})

    for r in rows:
        r['snr'] = np.sqrt(r['eff_hrs']) * r['weight']
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
    rows.sort(key=(lambda r: -r['snr']) if opts.sort == 'snr' else (lambda r: (not r['bright'], -r['snr'])))
    if opts.top is not None:
        rows = rows[:opts.top]
    return rows, n_all


CSV_HEADER = ['rank', 'hst_date', 'target', 'moon', 'illum', 'ut_mid', 't14_hrs', 'baseline_pre_hrs',
              'baseline_post_hrs', 'ut_start', 'ut_end', 'window_hrs', 'airmass_transit_mean',
              'airmass_transit_max', 'airmass_max', 'rel_snr']


def table_lines(rows, first):
    lines = [first,
             f"{'rank':>4} {'Date':10} {'target':14} {'moon':6} {'illum':>5} {'UT mid':>6} {'T14':>5} "
             f"{'pre':>5} {'post':>5} {'UT start':>8} {'UT end':>7} {'hrs':>5} {'X_tr':>5} {'Xmax':>5} {'S/N':>5}"]
    for i, r in enumerate(rows):
        lines.append(
            f"{i + 1:4d} {r['date']:10} {r['name'][:14]:14} {'bright' if r['bright'] else 'dark':6} "
            f"{r['illum']:5.2f} {r['tc']:6.2f} {r['t14']:5.2f} {r['pre']:5.2f} {r['post']:5.2f} "
            f"{r['ut0']:8.2f} {r['ut1']:7.2f} {r['hrs']:5.2f} {r['X_mean']:5.2f} {r['X_max']:5.2f} {r['snr']:5.2f}")
    return lines


def write_transits(path, rows, first, command=None):
    """Ranked transit windows to `path` (CSV) and the table, with the command, to .txt."""
    command = command if command is not None else ' '.join(sys.argv)
    with open(os.path.splitext(path)[0] + '.txt', 'w') as fh:
        fh.write(command + '\n' + '\n'.join(table_lines(rows, first)) + '\n')
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(CSV_HEADER)
        for i, r in enumerate(rows):
            w.writerow([i + 1, r['date'], r['name'], 'bright' if r['bright'] else 'dark', f"{r['illum']:.2f}",
                        f"{r['tc']:.3f}", f"{r['t14']:.3f}", f"{r['pre']:.2f}", f"{r['post']:.2f}",
                        f"{r['ut0']:.2f}", f"{r['ut1']:.2f}", f"{r['hrs']:.2f}", f"{r['X_mean']:.3f}",
                        f"{r['X_in_max']:.3f}", f"{r['X_max']:.3f}", f"{r['snr']:.3f}"])
    return [path, os.path.splitext(path)[0] + '.txt']
