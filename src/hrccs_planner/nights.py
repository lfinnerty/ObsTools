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
from dataclasses import dataclass

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


def observable_targets(targets, site, date, opts=None, log=print):
    """(night, [(target, altaz, velocity, draw_mask)]) for the targets shown on `date`."""
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
        if not (np.abs(kpmax - kpmin) > opts.delta_v_min and decreasing):
            continue
        draw = np.isfinite(kps) & (alt > opts.plot_alt_min)
        if opts.require_conjunction and not (kps[draw][0] > 0 and kps[draw][-1] < 0):
            continue
        shown.append((t, altaz, kps, draw))
    return night, shown


def plot_night(night, shown, site, path):
    """Two-panel figure (airmass; planet velocity) for one night, saved to `path`."""
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    hours = night.hours
    sunalt = night.sun_alt
    sun_secz = night.sun_altaz.secz
    moon_secz = night.moon_altaz.secz
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16, 12))
    for t, altaz, kps, draw in shown:
        up = altaz.alt.deg > 0
        ax[0].plot(hours[up], altaz.secz[up], label=t.name)
        ax[1].plot(hours[draw], kps[draw], label=t.name)
        ax[1].text(hours[draw][-1], kps[draw][-1], t.name)
    ax[0].plot(hours, sun_secz, color='y', label='Sun')
    ax[0].plot(hours, moon_secz, color='k', linestyle=':', label='Moon')
    for a, (lo, hi) in zip(ax, [(0, 90), (-300, 300)]):
        a.fill_between(hours, lo, hi, sunalt < 0., color='0.8', zorder=0)
        a.fill_between(hours, lo, hi, sunalt < -18., color='0.6', zorder=0)
        a.xaxis.set_major_locator(ticker.FixedLocator(np.arange(hours[0], hours[-1] + 0.01, 1)))
        a.grid(visible=True)
        a.set_xlabel('UT time [hr]')
    ax[1].axhline(0, color='k', linestyle='--')
    ax[0].legend()
    sunset = np.min(hours[sunalt < 0.]) - 1
    sunrise = np.max(hours[sunalt < 0.]) + 1
    ax[0].set_xlim(sunset, sunrise)
    ax[1].set_xlim(sunset, sunrise)
    ax[0].set_ylim([1.0, 3.0])
    ax[1].set_ylim(-200, 200)
    ax[0].set_ylabel('Airmass')
    ax[1].set_ylabel(r'$v_{pl}$ [km/s]')
    ax[0].set_title(f'{site.name}, night of {night.date} (local date)')
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def plan_nights(targets, site, dates, outdir, opts=None, log=print):
    """Plot every night with at least one observable target; returns the saved paths."""
    saved = []
    for date in dates:
        night, shown = observable_targets(targets, site, date, opts, log)
        if shown:
            path = os.path.join(outdir, f'{date}_{site.label}.png')
            log(f'Making plot for {date} at {site.label}: {", ".join(t.name for t, *_ in shown)}')
            plot_night(night, shown, site, path)
            saved.append(path)
    return saved
