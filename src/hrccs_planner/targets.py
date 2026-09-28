"""Target lists: one CSV row per planet, read and written by column name.

Required columns: Name, RA, Dec, Period (day), T transit/inf conj (JD),
Kp max (km/s). The other standard columns are optional; any extra columns are
kept in ``Target.extra`` and written back unchanged. See docs/targetlist_format.md.
"""
import csv
import math
import os
from dataclasses import dataclass, field
from typing import Optional

from astropy import units as u
from astropy.coordinates import SkyCoord

# (attribute, standard header) in the standard column order; the first 14
# columns of every target list, as read by the legacy scripts by position
STANDARD_COLUMNS = [
    ('name', 'Name'),
    ('ra', 'RA'),
    ('dec', 'Dec'),
    ('period', 'Period (day)'),
    ('a_au', 'a (AU)'),
    ('t0', 'T transit/inf conj (JD)'),
    ('status', 'Obs status'),
    ('eccentric', 'Eccentric?'),
    ('kp_max', 'Kp max (km/s)'),
    ('transiting', 'Transiting?'),
    ('kmag', 'Kmag (<9)'),
    ('e', 'e'),
    ('omega', 'omega'),
    ('inc', 'inc'),
]
STANDARD_HEADER = [h for _, h in STANDARD_COLUMNS]
REQUIRED = ['name', 'ra', 'dec', 'period', 't0', 'kp_max']

# other headers accepted for the standard columns (compared case-insensitively)
ALIASES = {
    'name': ['target', 'planet'],
    'period': ['period', 'period (d)', 'p (d)'],
    'a_au': ['a', 'a (au)'],
    't0': ['t0', 't0 (jd)', 'tc', 'tc (jd)', 't transit (jd)'],
    'status': ['kpic obs status', 'status'],
    'eccentric': ['eccentric'],
    'kp_max': ['kp max', 'kp (km/s)', 'kp'],
    'transiting': ['transiting'],
    'kmag': ['kmag', 'k'],
    'omega': ['omega (deg)', 'w'],
    'inc': ['inclination', 'i', 'inc (deg)'],
}

TRUE, FALSE = {'y', 'yes', 'true', 't', '1'}, {'n', 'no', 'false', 'f', '0', 'none', ''}


class TargetListError(ValueError):
    pass


@dataclass
class Target:
    name: str
    ra: str                  # sexagesimal hours (hh:mm:ss) or degrees
    dec: str                 # sexagesimal degrees or degrees
    period: float            # days
    t0: float                # transit / inferior conjunction (JD or BJD)
    kp_max: float            # total orbital velocity 2 pi a / P [km/s]
    a_au: Optional[float] = None
    status: str = ''
    eccentric: bool = False
    transiting: Optional[bool] = None
    kmag: float = math.nan
    e: float = 0.
    omega: float = 0.        # argument of periastron [deg]
    inc: Optional[float] = None
    extra: dict = field(default_factory=dict)

    @property
    def coord(self):
        unit = (u.deg, u.deg) if _is_number(self.ra) else (u.hourangle, u.deg)
        return SkyCoord(ra=self.ra, dec=self.dec, unit=unit)

    @property
    def host(self):
        """Host-star name: the planet name without a trailing planet letter."""
        parts = self.name.split()
        if len(parts) > 1 and len(parts[-1]) == 1 and parts[-1].islower():
            return ' '.join(parts[:-1])
        return self.name


def _is_number(text):
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


def _float(value, default=None):
    if value is None or str(value).strip() == '':
        return default
    return float(value)


def _bool(value, default=None):
    text = str(value).strip().lower() if value is not None else ''
    if text in TRUE:
        return True
    if text in FALSE:
        return False if text or default is None else default
    raise ValueError(f"expected Y/N, got '{value}'")


def _column_map(header):
    """{attribute: header} for the standard columns present in `header`."""
    lower = {h.strip().lower(): h for h in header if h is not None}
    out = {}
    for attr, std in STANDARD_COLUMNS:
        for cand in [std.lower()] + ALIASES.get(attr, []):
            if cand in lower:
                out[attr] = lower[cand]
                break
    return out


def read_targets(path):
    """List of Target from a target-list CSV (columns matched by name)."""
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames or []
        cols = _column_map(header)
        missing = [dict(STANDARD_COLUMNS)[a] for a in REQUIRED if a not in cols]
        if missing:
            raise TargetListError(f"{path}: missing column(s) {', '.join(missing)} "
                                  f"(header: {', '.join(header)})")
        standard = set(cols.values())
        targets = []
        for line, row in enumerate(reader, start=2):
            if not any((v or '').strip() for v in row.values() if isinstance(v, str)):
                continue
            get = lambda a, row=row: row.get(cols[a]) if a in cols else None
            try:
                e = _float(get('e'), 0.)
                ecc = get('eccentric')
                t = Target(
                    name=get('name').strip(), ra=get('ra').strip(), dec=get('dec').strip(),
                    period=_float(get('period')), t0=_float(get('t0')), kp_max=_float(get('kp_max')),
                    a_au=_float(get('a_au')), status=(get('status') or '').strip(),
                    eccentric=_bool(ecc) if ecc is not None and ecc.strip() else e > 0,
                    transiting=_bool(get('transiting'), None) if get('transiting') else None,
                    kmag=_float(get('kmag'), math.nan), e=e, omega=_float(get('omega'), 0.),
                    inc=_float(get('inc')),
                    extra={k: v for k, v in row.items() if k and k not in standard})
            except (TypeError, ValueError) as err:
                raise TargetListError(f'{path}, line {line} ({row.get(cols["name"])}): {err}') from None
            if t.period is None or t.t0 is None or t.kp_max is None:
                raise TargetListError(f'{path}, line {line} ({t.name}): Period, T0 and Kp max are required')
            targets.append(t)
    return targets


def _fmt(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ''
    if isinstance(value, bool):
        return 'Y' if value else 'N'
    return str(value)


def write_targets(path, targets, extra_columns=None):
    """Write targets with the standard columns followed by `extra_columns`
    (default: every extra column present in any target, in first-seen order)."""
    if extra_columns is None:
        extra_columns = list(dict.fromkeys(k for t in targets for k in t.extra))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(STANDARD_HEADER + list(extra_columns))
        for t in targets:
            row = [_fmt(getattr(t, a)) for a, _ in STANDARD_COLUMNS]
            row[STANDARD_HEADER.index('Eccentric?')] = 'Yes' if t.eccentric else 'No'
            w.writerow(row + [t.extra.get(c, '') for c in extra_columns])
