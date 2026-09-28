"""Sun, Moon and target positions over one night."""
from dataclasses import dataclass

import numpy as np
from astropy.coordinates import AltAz, get_body, get_sun


@dataclass
class Night:
    """Sun and Moon over the time grid of one night (see sites.Site.night_times)."""
    date: str
    times: object            # astropy Time array
    hours: np.ndarray        # UT hours relative to date 23:59:59 UTC
    frame: AltAz
    sun_altaz: object        # SkyCoord in `frame`
    moon_altaz: object       # SkyCoord in `frame`

    @property
    def sun_alt(self):
        return self._sun_alt

    def __post_init__(self):
        self._sun_alt = self.sun_altaz.alt.deg

    def dark(self, twilight=-18.):
        """Mask of times with the Sun below `twilight` degrees."""
        return self.sun_alt < twilight

    def moon_illumination(self, twilight=-18.):
        """Illuminated fraction of the Moon in the middle of the night."""
        night = self.dark(twilight)
        tmid = self.times[night][len(self.times[night]) // 2] if np.any(night) else self.times[len(self.times) // 2]
        elong = get_body('moon', tmid).separation(get_sun(tmid)).rad
        return (1 - np.cos(elong)) / 2.


def night_sky(site, date, n=1441):
    """Night for `site` starting on local date `date`, sampled at n points over 24 h
    (n = 1441 gives 1-minute steps)."""
    times, hours = site.night_times(date, n)
    frame = AltAz(obstime=times, location=site.location)
    sun_altaz = get_sun(times).transform_to(frame)
    moon_altaz = get_body('moon', times).transform_to(frame)
    return Night(date, times, hours, frame, sun_altaz, moon_altaz)


def airmass(alt_deg):
    """Plane-parallel airmass 1/sin(alt), with alt clipped to [1, 90] deg."""
    return 1. / np.sin(np.radians(np.clip(alt_deg, 1., 90.)))


def moon_too_close(target_altaz, night, min_sep_deg):
    """True if the Moon comes within min_sep_deg of the target while the target
    is above the horizon between sunset and sunrise."""
    up = (night.sun_alt < 0.) & (target_altaz.alt.deg > 0.)
    return bool(np.any(target_altaz.separation(night.moon_altaz).deg[up] < min_sep_deg))


def dates_between(start, end):
    """ISO dates from start through end, inclusive."""
    from datetime import date, timedelta
    d0, d1 = date.fromisoformat(start), date.fromisoformat(end)
    if d0 > d1:
        raise ValueError('start date must not be after end date')
    return [(d0 + timedelta(days=i)).isoformat() for i in range((d1 - d0).days + 1)]
