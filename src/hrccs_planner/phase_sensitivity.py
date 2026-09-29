"""How much ephemeris (phase) uncertainty costs dayside-emission HRCCS windows.

A window is planned from the predicted ephemeris and observed at fixed clock
times. If the true conjunction is off by delta (in phase), the observation still
happens at those times, but it covers phases shifted by delta. When the HRCCS
analysis refits the conjunction, delta changes only what is covered: the phase
curve f and the planet-velocity range, so the S/N of `windows` becomes

    S/N(delta) ~ sqrt(sum(f(phase + delta)^2 w dt) * delta_v(phase + delta)).

For each target's best windows this module averages S/N(delta) / S/N(0) over
delta ~ N(0, sigma_phase) and reports the mean and 10th percentile, for one or
more uncertainty scenarios. For an eccentric orbit, delta shifts the time since
transit; the dayside, eclipse, phase curve and velocities then follow the orbit as
in `windows` (ephemeris.eccentric_orbit). See docs/method.md.
"""
import csv
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
from astropy.time import Time

from .ephemeris import (
    eccentric_orbit,
    eclipse_phase,
    lambert_phase,
    orbital_phase,
    planet_rv_circular,
    propagate_conjunction,
)
from .visibility import airmass as airmass_of
from .visibility import night_sky
from .windows import rank_windows

# Phase shifts over which S/N(delta) is evaluated [cycles]; sigma_phase beyond
# ~0.06 cycles is truncated at the grid edge.
DELTA_GRID = np.linspace(-0.25, 0.25, 1001)
LOW_PERCENTILE = 10


@dataclass
class Scenario:
    """One ephemeris-uncertainty scenario.

    ``column`` names a target-list column with sigma_t [h] at the observing epoch.
    Without it, sigma_t comes from the T0 and period uncertainty columns,
    propagated to each window (ephemeris.propagate_conjunction).
    """

    name: str
    column: Optional[str] = None


def sigma_phase(target, scenario, jd, t0_err_column='T err (d)', period_err_column='P err (d)'):
    """Phase uncertainty [cycles] of `target` at time `jd` under `scenario` (NaN if unknown)."""
    if scenario.column is not None:
        try:
            return float(target.extra[scenario.column]) / 24. / target.period
        except (KeyError, TypeError, ValueError):
            return np.nan
    try:
        t0_err = float(target.extra.get(t0_err_column) or 0.)
        period_err = float(target.extra.get(period_err_column) or 0.)
    except ValueError:
        return np.nan
    return propagate_conjunction(target.t0, t0_err, target.period, period_err, jd)[3]


def snr_versus_shift(target, row, site, opts, deltas=DELTA_GRID):
    """S/N(delta) / S/N(0) for a window from `rank_windows`, at its planned clock times."""
    night = night_sky(site, row['date'], opts.samples)
    altaz = target.coord.transform_to(night.frame)
    alt = altaz.alt.deg
    X = airmass_of(alt)
    wam = 10 ** (-0.4 * opts.airmass_k * (X - 1.)) * X ** (-opts.seeing_exponent)
    phase = orbital_phase(night.times.jd, target.t0, target.period)
    above = alt > opts.alt_min
    if site.min_altitude is not None:
        above &= alt > site.min_altitude(altaz.az.deg)
    hours = night.hours
    ecl_halfwidth = 0.5 * opts.duration / 24. / target.period
    if target.eccentric:
        # as rank_windows: dayside by orbital phase, the eclipse around its own time
        geo = eccentric_orbit(phase, target.e, target.omega)[0] % 1
        out_of_eclipse = np.abs((phase - eclipse_phase(target.e, target.omega) + 0.5) % 1 - 0.5) >= ecl_halfwidth
    else:
        geo = phase
        out_of_eclipse = np.abs(phase - 0.5) >= ecl_halfwidth
    # The same samples rank_windows used: usable time within the window's UT span.
    mask = (night.dark(opts.twilight) & above & (hours >= row['ut0'] - 1e-9) & (hours <= row['ut1'] + 1e-9)
            & (geo > 0.25) & (geo < 0.75) & out_of_eclipse)
    shifted = (phase[mask][np.newaxis, :] + np.asarray(deltas)[:, np.newaxis]) % 1
    kp = target.kp_max * np.sin(np.radians(opts.inclination))
    if target.eccentric:
        # delta shifts the time since transit; the orbit gives the geometry and velocity
        uphase, _, rv = eccentric_orbit(shifted, target.e, target.omega)
        signal = np.sum(lambert_phase(uphase % 1, opts.inclination) ** 2 * wam[mask], axis=1)
        velocity = kp * rv
    else:
        signal = np.sum(lambert_phase(shifted, opts.inclination) ** 2 * wam[mask], axis=1)
        velocity = planet_rv_circular(shifted, kp)
    snr = np.sqrt(signal * (velocity.max(axis=1) - velocity.min(axis=1)))
    return snr / snr[np.argmin(np.abs(deltas))]


def expected_ratio(ratio, sigma, deltas=DELTA_GRID, percentile=LOW_PERCENTILE):
    """Mean and low percentile of S/N(delta)/S/N(0) for delta ~ N(0, sigma) [cycles]."""
    if not np.isfinite(sigma):
        return np.nan, np.nan
    if sigma <= 0:
        return 1., 1.
    weights = np.exp(-0.5 * (np.asarray(deltas) / sigma) ** 2)
    weights /= weights.sum()
    order = np.argsort(ratio)
    cdf = np.cumsum(weights[order])
    low = ratio[order][min(np.searchsorted(cdf, percentile / 100.), len(order) - 1)]
    return float(np.sum(weights * ratio)), float(low)


def phase_sensitivity(targets, site, dates, opts, scenarios, top=5, t0_err_column='T err (d)',
                      period_err_column='P err (d)', log=print):
    """Expected S/N loss from phase uncertainty for each target's best windows.

    Returns (summary, windows): one dict per target and scenario, averaging its
    ``top`` best windows (by relative S/N), and one dict per window and scenario.
    """
    rows, _ = rank_windows(targets, site, dates, opts, log=log)
    by_name = {t.name: t for t in targets}
    summary, windows = [], []
    for name in dict.fromkeys(r['name'] for r in rows):
        target = by_name[name]
        best = sorted((r for r in rows if r['name'] == name), key=lambda r: -r['snr'])[:top]
        ratios = [snr_versus_shift(target, row, site, opts) for row in best]
        for scenario in scenarios:
            per_window = []
            for row, ratio in zip(best, ratios):
                jd = Time(row['date'] + 'T12:00:00', scale='utc').jd + 0.5 + row['ut0'] / 24.
                sig = sigma_phase(target, scenario, jd, t0_err_column, period_err_column)
                mean, low = expected_ratio(ratio, sig)
                per_window.append((row, sig, mean, low))
                windows.append({'name': name, 'scenario': scenario.name, 'date': row['date'],
                                'ut0': row['ut0'], 'ut1': row['ut1'], 'ph0': row['ph0'], 'ph1': row['ph1'],
                                'rel_snr': row['snr'], 'sigma_t_h': sig * target.period * 24.,
                                'mean_ratio': mean, 'low_ratio': low})
            summary.append({
                'name': name, 'scenario': scenario.name, 'period': target.period, 'n_windows': len(best),
                'sigma_t_h': float(np.nanmedian([w[1] for w in per_window])) * target.period * 24.,
                'mean_ratio': float(np.nanmean([w[2] for w in per_window])),
                'low_ratio': float(np.nanmean([w[3] for w in per_window])),
            })
    return summary, windows


def table_lines(summary, scenarios, percentile=LOW_PERCENTILE):
    """Printable table: one line per target, with each scenario's S/N ratios."""
    names = [s.name for s in scenarios]
    group = 24  # width of one scenario's column group
    top = f"{'':32s}" + ''.join(f' | {name[:group - 1]:^{group - 1}s}' for name in names)
    head = f"{'target':16s} {'P [d]':>7s} {'windows':>7s}"
    head += ''.join(f" | {'sig_t [h]':>9s} {'<S/N>':>6s} {f'p{percentile}':>6s}" for _ in names)
    if len(names) > 1:
        top += f" | {'gain':>6s}"
        head += f" | {'<S/N>':>6s}"
    lines = [top, head]
    for target in dict.fromkeys(s['name'] for s in summary):
        rows = {s['scenario']: s for s in summary if s['name'] == target}
        first = rows[names[0]]
        line = f"{target[:16]:16s} {first['period']:7.3f} {first['n_windows']:7d}"
        for name in names:
            r = rows[name]
            line += f" | {r['sigma_t_h']:9.2f} {r['mean_ratio']:6.3f} {r['low_ratio']:6.3f}"
        if len(names) > 1:
            gains = [100 * (rows[n]['mean_ratio'] / first['mean_ratio'] - 1) for n in names[1:]]
            line += ' | ' + ' '.join(f'{g:+5.1f}%' for g in gains)
        lines.append(line)
    if len(names) > 1:
        lines.append(f"(gain: change in mean S/N relative to '{names[0]}')")
    return lines


def write_summary(path, summary, windows):
    """Write the per-target summary to `path` and the per-window rows next to it."""
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['target', 'scenario', 'period_d', 'n_windows', 'sigma_t_h', 'mean_snr_ratio', 'p10_snr_ratio'])
        for s in summary:
            writer.writerow([s['name'], s['scenario'], f"{s['period']:.6f}", s['n_windows'],
                             f"{s['sigma_t_h']:.3f}", f"{s['mean_ratio']:.4f}", f"{s['low_ratio']:.4f}"])
    root, ext = os.path.splitext(path)
    window_path = f'{root}_windows{ext or ".csv"}'
    with open(window_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['target', 'scenario', 'hst_date', 'ut_start', 'ut_end', 'phase_start', 'phase_end',
                         'rel_snr', 'sigma_t_h', 'mean_snr_ratio', 'p10_snr_ratio'])
        for w in windows:
            writer.writerow([w['name'], w['scenario'], w['date'], f"{w['ut0']:.2f}", f"{w['ut1']:.2f}",
                             f"{w['ph0']:.3f}", f"{w['ph1']:.3f}", f"{w['rel_snr']:.3f}", f"{w['sigma_t_h']:.3f}",
                             f"{w['mean_ratio']:.4f}", f"{w['low_ratio']:.4f}"])
    return [path, window_path]
