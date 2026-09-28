"""Telescope starlists from SIMBAD positions and proper motions.

SIMBAD gives the proper motion in RA as mu_alpha* = mu_alpha cos(dec) [mas/yr].
Both formats here take the rate of change of the RA coordinate in seconds of
time per year, mu_alpha* / (15000 cos dec), and mu_dec in arcsec/yr.
"""
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from .catalogs import query_simbad

KECK_COLWIDTH = 16
MAGELLAN_HEADER = ('# RA Dec equinox RApm Decpm offset rot RA_probe1 Dec_probe1 equinox RA_probe2 '
                   'Dec_probe2 equinox pm_epoch\n###name hh:mm:ss.s sdd:mm:ss yyyy.0 s.ss s.ss angle mode '
                   'hh:mm:ss.s sdd:mm:ss yyyy.0 hh:mm:ss.s sdd:mm:ss yyyy.0 yyyy.0\n')


def _mag(entry, band):
    v = entry[band] if band in entry.colnames else None
    return None if v is None or np.ma.is_masked(v) else float(v)


def _pm(entry, coord):
    """(RA rate [s of time/yr], Dec rate [arcsec/yr])."""
    pmra = 0. if np.ma.is_masked(entry['pmra']) else float(entry['pmra'])
    pmdec = 0. if np.ma.is_masked(entry['pmdec']) else float(entry['pmdec'])
    return pmra / 1e3 / (15. * np.cos(coord.dec.radian)), pmdec / 1e3


def keck_lines(names, rows):
    """Keck starlist lines: name, RA, Dec, equinox, pmra=, pmdec=, and R (or V) magnitude."""
    lines = []
    for name, entry in zip(names, rows):
        c = SkyCoord(ra=float(entry['ra']) * u.deg, dec=float(entry['dec']) * u.deg)
        ra = c.ra.to_string(unit=u.hourangle, sep=' ', precision=3, pad=True)[:12]
        dec = c.dec.to_string(unit=u.deg, sep=' ', precision=2, pad=True, alwayssign=True)[:12]
        pmra, pmdec = _pm(entry, c)
        s = name[:15].ljust(KECK_COLWIDTH) + ra.ljust(KECK_COLWIDTH) + dec.ljust(KECK_COLWIDTH) + '2000 '
        s += f'pmra={np.round(pmra, 4)}'.ljust(KECK_COLWIDTH) + f'pmdec={np.round(pmdec, 4)}'.ljust(KECK_COLWIDTH)
        r, v = _mag(entry, 'R'), _mag(entry, 'V')
        if r is not None:
            s += f'Rmag={np.round(r, 2)}'.ljust(KECK_COLWIDTH)
        elif v is not None:
            s += f'Vmag={np.round(v, 2)}'.ljust(KECK_COLWIDTH)
        lines.append(s.rstrip())
    return lines


def magellan_lines(names, rows, rot_angle='0.0', rot_mode='GRV', epoch='2017.5'):
    """Magellan (TCS catalog format) lines, numbered from 001, with the V magnitude as a comment."""
    lines = []
    for i, (name, entry) in enumerate(zip(names, rows)):
        c = SkyCoord(ra=float(entry['ra']) * u.deg, dec=float(entry['dec']) * u.deg)
        ra = c.ra.to_string(unit=u.hourangle, sep=':', precision=2, pad=True)
        dec = c.dec.to_string(unit=u.deg, sep=':', precision=2, pad=True, alwayssign=True)
        pmra, pmdec = _pm(entry, c)
        fields = [str(i + 1).zfill(3), name[:15].replace(' ', '_'), ra, dec, '2000.0', str(np.round(pmra, 3)),
                  str(np.round(pmdec, 3)), rot_angle, rot_mode, '0', '0', '2000.0', '0', '0', '2000.0', epoch]
        s = ' '.join(fields)
        v = _mag(entry, 'V')
        if v is not None:
            s += f' # V={np.round(v, 1)}'
        lines.append(s)
    return lines


def write_starlist(path, names, fmt='keck', workspace=None, refresh=False, log=print):
    """Query SIMBAD for `names` and write a starlist in `fmt` ('keck' or 'magellan')."""
    rows = query_simbad(names, workspace=workspace, refresh=refresh, log=log)
    if fmt == 'keck':
        text = '\n'.join(keck_lines(names, rows)) + '\n'
    elif fmt == 'magellan':
        text = MAGELLAN_HEADER + '\n'.join(magellan_lines(names, rows)) + '\n'
    else:
        raise ValueError(f"unknown starlist format '{fmt}' (keck, magellan)")
    with open(path, 'w') as fh:
        fh.write(text)
    return text
