"""Nightly planning plots: airmass and planet velocity of every target that is
observable on a night.

A target is shown on a night if (all thresholds are options):
    - the Moon stays at least moon_sep_min from it while it is up at night;
    - it is above alt_limit (30 deg) for more than min_samples time samples
      between sunset and sunrise (Sun below 0 deg);
    - its planet velocity changes by more than delta_v_min while above alt_limit;
    - the velocity decreases over the night (the planet is on the dayside half of
      the orbit, around secondary eclipse, for circular orbits);
    - with require_conjunction, the velocity crosses zero (from + to -).
Eccentric orbits use the Keplerian velocity curve.
"""
import os
from dataclasses import dataclass, replace

import numpy as np

from .ephemeris import planet_rv_eccentric
from .visibility import moon_too_close, night_sky


@dataclass
class NightOptions:
    alt_limit: float = 30.
    min_samples: int = 4
    delta_v_min: float = 30.
    plot_alt_min: float = 15.       # velocity curves are drawn where the target is above this
    require_conjunction: bool = False
    moon_sep_min: float = 10.
    inclination: float = 90.
    samples: int = 100


def planet_velocity(target, jd, kp):
    """Planet velocity [km/s] over one night (legacy phase convention: the
    number of orbits since T0 at the start of the night is removed)."""
    phase = (jd - target.t0) / target.period
    phase -= int(phase[0])
    if not target.eccentric:
        return kp * np.sin(phase * 2 * np.pi)
    return planet_rv_eccentric(phase, kp, target.e, target.omega)


def observable_targets(targets, site, date, opts=None, log=print, always=()):
    """(night, [(target, altaz, velocity, draw_mask)]) for the targets shown on `date`.
    Targets named in `always` (e.g. with a transit window to shade) skip the
    velocity-change, dayside and conjunction criteria."""
    opts = opts or NightOptions()
    night = night_sky(site, date, opts.samples)
    sunalt = night.sun_alt
    shown = []
    for t in targets:
        kp = t.kp_max * np.sin(np.radians(opts.inclination))
        kps = planet_velocity(t, night.times.jd, kp)
        altaz = t.coord.transform_to(night.frame)
        alt = altaz.alt.deg
        if moon_too_close(altaz, night, opts.moon_sep_min):
            up = (sunalt < 0.) & (alt > 0.)
            sep = altaz.separation(night.moon_altaz).deg[up]
            log(f'Skipping {t.name} on {date} - min moon separation {np.min(sep):.1f} deg')
            continue
        if np.sum(alt[sunalt < 0.] > opts.alt_limit) <= opts.min_samples:
            continue
        kps = np.where(sunalt > 0., np.nan, kps)
        kpmin = np.nanmin(kps[alt > opts.alt_limit])
        kpmax = np.nanmax(kps[alt > opts.alt_limit])
        kpclean = kps[~np.isnan(kps)]
        decreasing = kpclean[-1] < kpclean[0]
        draw = np.isfinite(kps) & (alt > opts.plot_alt_min)
        if t.name not in always:
            if not (np.abs(kpmax - kpmin) > opts.delta_v_min and decreasing):
                continue
            if opts.require_conjunction and not (kps[draw][0] > 0 and kps[draw][-1] < 0):
                continue
        shown.append((t, altaz, kps, draw))
    return night, shown


def host_key(target):
    """Targets orbiting the same star (same position to ~1 arcsec) share a key."""
    c = target.coord
    return (round(c.ra.deg, 3), round(c.dec.deg, 3))


def star_label(host):
    """Primary-star label: 'TOI-1408' -> 'TOI-1408A', 'HD 189733' -> 'HD 189733A', '55 Cnc' -> '55 Cnc A'."""
    return host + ('A' if host[-1:].isdigit() else ' A')


def read_shading(paths):
    """Observing windows to shade, from ranked-window or transit-window CSVs
    (`hrccs-plan windows/transits --csv`): {(date, target): [(ut0, ut1, t1, t4)]},
    with the transit contacts t1, t4 (UT h) for transit windows, else None."""
    import csv
    out = {}
    for path in paths:
        with open(path, newline='') as fh:
            for r in csv.DictReader(fh):
                t1 = t4 = None
                if r.get('ut_mid') and r.get('t14_hrs'):
                    t1 = float(r['ut_mid']) - float(r['t14_hrs']) / 2.
                    t4 = float(r['ut_mid']) + float(r['t14_hrs']) / 2.
                out.setdefault((r['hst_date'], r['target']), []).append(
                    (float(r['ut_start']), float(r['ut_end']), t1, t4))
    return out


def parse_window(text):
    """'HH:MM-HH:MM' or 'H.H-H.H' (UT on the night's grid; hours < 12 are after
    00 UT, i.e. relative to date 24:00 UT) -> (ut0, ut1) in hours."""
    def h(x):
        x = x.strip()
        if ':' in x:
            hh, mm = x.split(':')
            return int(hh) + int(mm) / 60.
        return float(x)
    a, b = text.split('-') if text.count('-') == 1 else text.rsplit('-', 1)
    ut0, ut1 = h(a), h(b)
    if ut1 < ut0:
        ut0 -= 24.
    return ut0, ut1


def plot_night(night, shown, site, path, windows=None, compact=False, title=None, site_label=None):
    """Two-panel figure (airmass; planet velocity) for one night, saved to `path`.

    windows -- {target name: [(ut0, ut1, t1, t4)]} observing windows to shade (t1, t4:
               transit contacts, shaded darker, or None)
    compact -- a small figure for proposals and papers: shared UT axis, the dark part
               of the night only, the airmass-2 limit marked, the Moon drawn only
               above the horizon and the planet velocity highlighted in the window
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    windows = windows or {}
    hours = night.hours
    sunalt = night.sun_alt
    moon_alt = night.moon_altaz.alt.deg
    moon_secz = np.where(moon_alt > 0, night.moon_altaz.secz, np.nan)
    if compact:
        plt.rcParams.update({'font.size': 10})
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6.4, 5.2), sharex=True,
                               gridspec_kw={'hspace': 0.08})
    else:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16, 12))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    shaded = False
    # one airmass curve per host star: planets of the same star share it, labelled
    # with the host name (+ 'A') when more than one of its planets is shown
    hosts = {}
    for t, altaz, *_ in shown:
        hosts.setdefault(host_key(t), []).append((t, altaz))
    for members in hosts.values():
        t, altaz = members[0]
        multi = len(members) > 1
        c = '0.25' if multi else colors[[x.name for x, *_ in shown].index(t.name) % len(colors)]
        up = altaz.alt.deg > 0
        ax[0].plot(hours[up], altaz.secz[up], color=c, label=star_label(t.host) if multi else t.name,
                   lw=1.6 if compact else 1.5)
    vel_handles = []
    for k, (t, _, kps, draw) in enumerate(shown):
        c = colors[k % len(colors)]
        line, = ax[1].plot(hours[draw], kps[draw], color=c, label=t.name, lw=1.2 if compact else 1.5,
                           alpha=0.6 if compact and t.name in windows else 1.)
        vel_handles.append(line)
        if not compact:
            ax[1].text(hours[draw][-1], kps[draw][-1], t.name)
        for ut0, ut1, t1, t4 in windows.get(t.name, []):
            lab = 'Observing window' if not shaded else None
            for a in ax:
                a.axvspan(ut0, ut1, facecolor=c, alpha=0.15, lw=0, zorder=0.5, label=lab if a is ax[0] else None)
                for x in (ut0, ut1):
                    a.axvline(x, color=c, lw=0.9, alpha=0.8, zorder=0.7)
                if t1 is not None:
                    a.axvspan(t1, t4, facecolor=c, alpha=0.25, lw=0, zorder=0.6,
                              label=('In transit' if not shaded else None) if a is ax[0] else None)
            shaded = True
            inw = draw & (hours >= ut0) & (hours <= ut1)
            ax[1].plot(hours[inw], kps[inw], color=c, lw=3 if compact else 2.5, solid_capstyle='butt')
    if not compact:
        ax[0].plot(hours, night.sun_altaz.secz, color='y', label='Sun')
    ax[0].plot(hours, moon_secz, color='k', linestyle=':', label='Moon')
    for a, (lo, hi) in zip(ax, [(0, 90), (-300, 300)]):
        a.fill_between(hours, lo, hi, sunalt < 0., color='0.94' if compact else '0.8', zorder=0, lw=0)
        a.fill_between(hours, lo, hi, sunalt < -18., color='0.86' if compact else '0.6', zorder=0, lw=0)
        a.xaxis.set_major_locator(ticker.FixedLocator(np.arange(np.ceil(hours[0]), hours[-1] + 0.01, 1)))
        a.grid(visible=True, alpha=0.4 if compact else 1.)
    ax[1].axhline(0, color='k', linestyle='--', lw=0.8 if compact else 1.5)
    if compact:
        dark = hours[sunalt < 0.]
        xlo, xhi = np.min(dark) - 0.2, np.max(dark) + 0.2
        ax[0].axhline(2., color='0.3', lw=0.8, ls='--')
        ax[0].set_ylim(2.6, 1.0)          # airmass increases downward
        ax[0].set_ylabel('Airmass')
        ax[1].set_xlabel('UT [h]')
        vmax = max(np.nanmax(np.abs(kps[draw])) for _, _, kps, draw in shown)
        vlim = 50. * np.ceil(1.15 * vmax / 50.)
        ax[1].set_ylim(-vlim, vlim)
        ax[1].set_ylabel(r'$v_\mathrm{pl}$ [km s$^{-1}$]')
        ax[0].legend(fontsize=8, loc='lower center', ncol=3, frameon=True, framealpha=0.9)
        if len(shown) > 1:
            ax[1].legend(handles=vel_handles, fontsize=8, loc='lower left', frameon=True, framealpha=0.9)
        ax[1].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x % 24:.0f}'))
        for a in ax:
            a.tick_params(direction='in', top=True, right=True)
        parts = []
        for members in hosts.values():
            letters = [t.name[len(t.host):].strip() for t, _ in members]
            if len(members) > 1 and all(len(x) == 1 for x in letters):
                parts.append(f"{members[0][0].host} {', '.join(letters[:-1])} and {letters[-1]}")
            else:
                parts.extend(t.name for t, _ in members)
        names = ', '.join(parts)
        ax[0].set_title(title if title is not None else
                        f'{names}: {site_label or site.name}, night of {night.date}', fontsize=11)
    else:
        for a in ax:
            a.set_xlabel('UT time [hr]')
        ax[0].legend()
        xlo, xhi = np.min(hours[sunalt < 0.]) - 1, np.max(hours[sunalt < 0.]) + 1
        ax[0].set_ylim([1.0, 3.0])
        ax[1].set_ylim(-200, 200)
        ax[0].set_ylabel('Airmass')
        ax[1].set_ylabel(r'$v_{pl}$ [km/s]')
        ax[0].set_title(title if title is not None else f'{site.name}, night of {night.date} (local date)')
    ax[0].set_xlim(xlo, xhi)
    ax[1].set_xlim(xlo, xhi)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, bbox_inches='tight', dpi=300 if compact else 100)
    plt.close(fig)
    if compact:
        plt.rcParams.update(plt.rcParamsDefault)


def plan_nights(targets, site, dates, outdir, opts=None, log=print, shading=None, compact=False,
                fmt='png', title=None, site_label=None, show_all=False):
    """Plot every night with at least one observable target; returns the saved paths.

    shading -- {(date, target): [(ut0, ut1, t1, t4)]} windows to shade (read_shading),
               or {None: [...]} to shade the same windows on every night
    """
    saved = []
    shading = shading or {}
    for date in dates:
        always = {name for (d, name) in shading if d == date} | (
            {t.name for t in targets} if None in shading or show_all else set())
        night, shown = observable_targets(targets, site, date, opts, log, always)
        if shown and compact:
            # same targets, redrawn on a 2-minute grid for smooth curves and exact twilight
            o = opts or NightOptions()
            names = {t.name for t, *_ in shown}
            fine = replace(o, samples=721, min_samples=int(round(o.min_samples * 7.2)))
            night, shown2 = observable_targets([t for t in targets if t.name in names], site, date, fine,
                                               lambda msg: None, always)
            shown = shown2 or shown
        if shown:
            wins = {t.name: shading.get((date, t.name), shading.get(None, [])) for t, *_ in shown}
            wins = {k: v for k, v in wins.items() if v}
            suffix = '_compact' if compact else ''
            path = os.path.join(outdir, f'{date}_{site.label}{suffix}.{fmt}')
            log(f'Making plot for {date} at {site.label}: {", ".join(t.name for t, *_ in shown)}'
                + (f' (shaded: {", ".join(wins)})' if wins else ''))
            plot_night(night, shown, site, path, wins, compact, title, site_label)
            saved.append(path)
    return saved
