"""Observatory sites and the time grid used for one night.

Each night is labelled by the local calendar date on which it starts (at
Maunakea this is the HST date). The time grid for that night is a 24-hour span
of UT times, starting at ``date 23:59:59 UTC + hours_start`` hours, so that it
covers the whole local night. ``hours_start`` is -3 for Maunakea and -6 for the
Chilean sites; for other sites it is derived from the longitude.

Built-in sites use fixed geocentric coordinates (those of the astropy site
registry), so no network access is needed. Any other astropy site name, e.g.
``"Apache Point"``, also works, via ``EarthLocation.of_site``.
"""
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time


@dataclass(frozen=True)
class Site:
    key: str                 # registry key, e.g. "maunakea"
    name: str                # human-readable name
    label: str               # short label used in plot file names, e.g. "MKO"
    xyz_m: tuple             # geocentric (x, y, z) [m]
    hours_start: float       # UT offset of the night grid start from date 23:59:59 UTC [h]
    min_altitude: Optional[Callable] = None  # optional pointing limit: alt_min(az_deg) [deg]

    @property
    def location(self):
        return EarthLocation.from_geocentric(*self.xyz_m, unit=u.m)

    def night_hours(self, n):
        """UT hours (relative to date 23:59:59 UTC) of an n-point grid over 24 h."""
        return np.linspace(self.hours_start, self.hours_start + 24., n)

    def night_times(self, date, n):
        """(Time array, hours) for the night starting on local date `date` (YYYY-MM-DD)."""
        hours = self.night_hours(n)
        return Time(date + 'T23:59:59', format='isot', scale='utc') + hours * u.hour, hours


def keck2_min_altitude(az_deg):
    """Keck II pointing limit [deg]: 36.8 deg over azimuths 185.3-332.8 deg (the
    Nasmyth deck), 18 deg elsewhere (https://www2.keck.hawaii.edu/inst/common/TelLimits.html)."""
    az = np.asarray(az_deg) % 360.
    return np.where((az > 185.3) & (az < 332.8), 36.8, 18.)


_KECK = (-5464487.817598869, -2492806.59108569, 2151240.1945184576)
_LCO = (1845655.4990534068, -5270856.2947176, -3075330.777606821)
_PARANAL = (1946404.3410388362, -5467644.290798524, -2642728.2014442487)
_GEMINI_S = (1820193.0684460273, -5208343.034275673, -3194842.5004834323)

SITES = {
    'maunakea': Site('maunakea', 'Maunakea (Keck)', 'MKO', _KECK, -3.),
    'keck2': Site('keck2', 'Keck II (with pointing limits)', 'MKO', _KECK, -3., keck2_min_altitude),
    'lco': Site('lco', 'Las Campanas Observatory', 'LCO', _LCO, -6.),
    'paranal': Site('paranal', 'Paranal', 'EELT', _PARANAL, -6.),
    'gemini-s': Site('gemini-s', 'Cerro Pachon (Gemini South)', 'Gem-S', _GEMINI_S, -6.),
}

# Telescope and instrument names accepted for each site (case-insensitive).
ALIASES = {
    'maunakea': ['keck', 'gemini-n', 'gemini_north', 'mko', 'irtf', 'subaru', 'cfht', 'kpic',
                 'igrins2', 'nirspec', 'hispec'],
    'keck2': [],
    'lco': ['magellan', 'las campanas', 'las campanas observatory', 'mike', 'winered'],
    'paranal': ['vlt', 'eelt', 'elt', 'eso', 'crires', 'crires+'],
    'gemini-s': ['gemini_south', 'soar', 'ctio', 'cerro pachon'],
}


def get_site(name):
    """Site for a registry key, an alias (telescope or instrument name) or an
    astropy site name."""
    key = name.strip().lower()
    for site_key, aliases in ALIASES.items():
        if key == site_key or key in aliases:
            return SITES[site_key]
    try:
        loc = EarthLocation.of_site(name)
    except Exception as err:
        known = ', '.join(sorted(SITES) + sorted(a for v in ALIASES.values() for a in v))
        raise ValueError(f"unknown site '{name}' (known: {known}; or any astropy site name)") from err
    # start the grid at (mean solar) local noon of the date, so it spans the night
    local_noon_ut = 12. - loc.lon.deg / 15.
    hours_start = float(np.floor(local_noon_ut - 24.))
    xyz = tuple(float(c.to_value(u.m)) for c in loc.to_geocentric())
    return Site(key, name, name.replace(' ', '_'), xyz, hours_start)
