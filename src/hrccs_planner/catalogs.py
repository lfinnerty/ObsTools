"""Batched catalog queries, cached on disk.

Results are saved in <workspace>/cache/ as ECSV files named by a hash of the
query, and reused instead of querying again. Failed or empty queries are not
cached. Pass refresh=True to query again (e.g. after a catalog update).
"""
import hashlib
import math
import os
import re

import numpy as np

from .targets import Target
from .workspace import subdir

AU_PER_DAY_TO_KM_S = 1731.45683633
EXOARCHIVE_COLUMNS = ('hostname,pl_name,ra,dec,pl_orbper,pl_orbsmax,pl_tranmid,'
                      'pl_orbeccen,pl_orblper,pl_orbincl,sy_kmag,tran_flag,pl_trandur')
T14_COLUMN = 'T14 (h)'  # transit duration from the archive (pl_trandur), used by `transits --duration-column`


def _cached(kind, key, query, workspace=None, refresh=False, log=print):
    from astropy.table import Table
    path = os.path.join(subdir('cache', workspace, create=True),
                        f'{kind}_{hashlib.sha1(key.encode()).hexdigest()[:16]}.ecsv')
    if os.path.exists(path) and not refresh:
        return Table.read(path, format='ascii.ecsv')
    log(f'Querying {kind} ...')
    table = query()
    if table is not None and len(table):
        Table(table).write(path, format='ascii.ecsv', overwrite=True)
    return table


def _value(row, column):
    value = row[column]
    if value is None or np.ma.is_masked(value):
        return None
    if hasattr(value, 'unit') and hasattr(value, 'value'):
        value = value.value
    try:
        value = value.item()
    except AttributeError:
        pass
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _name_variants(name):
    variants = [name]
    spaced = re.sub(r'(?<=\d)([A-Z])(?= [a-z]$)', r' \1', name)
    if spaced != name:
        variants.append(spaced)
    return variants


def query_exoarchive(system_names, workspace=None, refresh=False, log=print):
    """NASA Exoplanet Archive (pscomppars) rows for hosts or planets, one batched query."""
    def normalized(name):
        return re.sub(r'\s+', '', name).casefold()
    conditions = []
    for name in system_names:
        for v in _name_variants(name):
            e = v.replace("'", "''")
            conditions.append(f"hostname = '{e}' OR pl_name = '{e}'")
    where = ' OR '.join(conditions)

    def query():
        from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
        return NasaExoplanetArchive.query_criteria(table='pscomppars', select=EXOARCHIVE_COLUMNS, where=where)
    rows = _cached('exoarchive', where, query, workspace, refresh, log)
    found = {normalized(str(_value(r, c))) for r in rows for c in ('hostname', 'pl_name') if _value(r, c)}
    missing = [n for n in system_names if normalized(n) not in found]
    if missing:
        raise ValueError('No planets found in the NASA Exoplanet Archive for: ' + ', '.join(missing))
    return rows


def exoarchive_target(row, circular_below=0.):
    """Target from a NASA Exoplanet Archive pscomppars row. Kp max = 2 pi a / P
    (no eccentricity correction; the planners apply the Keplerian velocity curve).
    Orbits with e < circular_below are written as circular (e = 0, Eccentric? = No)."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c = SkyCoord(ra=float(_value(row, 'ra')) * u.deg, dec=float(_value(row, 'dec')) * u.deg)
    ra, dec = c.to_string('hmsdms', sep=':', precision=2, pad=True).split()
    per, a = _value(row, 'pl_orbper'), _value(row, 'pl_orbsmax')
    e = _value(row, 'pl_orbeccen')
    if e is not None and float(e) < circular_below:
        e = 0.
    kp = round(2 * math.pi * float(a) / float(per) * AU_PER_DAY_TO_KM_S) if per and a else None
    kmag = _value(row, 'sy_kmag')
    return Target(
        name=_value(row, 'pl_name'), ra=ra, dec=dec, period=per, t0=_value(row, 'pl_tranmid'), kp_max=kp,
        a_au=a, status='None', eccentric=e is not None and float(e) > 0,
        transiting=_value(row, 'tran_flag') == 1, kmag=kmag if kmag is not None else math.nan,
        e=e if e is not None else 0., omega=_value(row, 'pl_orblper'), inc=_value(row, 'pl_orbincl'),
        extra={T14_COLUMN: _value(row, 'pl_trandur') if _value(row, 'pl_trandur') is not None else ''})


def query_tic(tic_ids, workspace=None, refresh=False, log=print):
    """{TIC ID: row} from one batched TIC-8 query at MAST (Stassun et al. 2019)."""
    ids = sorted({int(t) for t in tic_ids})

    def query():
        from astroquery.mast import Catalogs
        return Catalogs.query_criteria(catalog='Tic', ID=ids)
    rows = _cached('tic', ','.join(map(str, ids)), query, workspace, refresh, log)
    return {int(r['ID']): r for r in rows}


def query_simbad(names, fields=('pmra', 'pmdec', 'V', 'R'), workspace=None, refresh=False, log=print):
    """SIMBAD rows (main_id, ra, dec in deg, and `fields`) for `names`, in order."""
    def query():
        from astroquery.simbad import Simbad
        s = Simbad()
        s.add_votable_fields(*fields)
        return s.query_objects(list(names))
    rows = _cached('simbad', '|'.join(names) + '#' + ','.join(fields), query, workspace, refresh, log)
    bad = [n for n, r in zip(names, rows) if not str(r['main_id']).strip()]
    if len(rows) != len(names) or bad:
        raise ValueError('SIMBAD found no object for: ' + ', '.join(bad or names))
    return rows
