"""
Make an observing target list (targetlists/<name>_targetlist.csv, the format
read by observing_planner.py and best_window.py) from the ranked candidates of
the TESS blind phase-curve search for non-transiting ultra-hot Jupiters
(~/Research/TESS: tess_blind_rank.py -> candidates_ranked*.csv).

Columns 0-13 follow the existing target lists (the planners index them by
position):
	Name, RA, Dec, Period (day), a (AU), T transit/inf conj (JD), Obs status,
	Eccentric?, Kp max (km/s), Transiting?, Kmag (<9), e, omega, inc
followed by extra columns from the search (tier, TESS amplitude, implied
radius, ephemeris errors, period aliases, notes).

- RA/Dec (J2000) and Kmag come from one batched TIC-8 query at MAST
  (Stassun et al. 2019).
- a is from Kepler's law with the stellar mass used in the ranking (TIC, or
  estimated from V, distance and Teff where the TIC has none); Kp max =
  2 pi a / P, i.e. sin i = 1, since the inclination of a non-transiting
  planet is unknown.
- Expected dayside contrast in the IGRINS H (1.49-1.80 um) and K (1.96-2.48
  um) bands: (Rp/R*)^2 B(T_day)/B(T_eff), photon-weighted blackbody band
  averages, with T_day from Cowan & Agol (2011) for Bond albedo 0.1 and
  redistribution 0.2 (the nominal model of the ranking) and Rp = --rp (default
  1.5 R_J). Thermal emission only; reflected light is negligible in the IR.
  best_window.py --weight-contrast uses the K-band column.
- T is the inferior conjunction (phase-curve minimum) projected to the first
  epoch after the ranking's follow-up date (default 2027-01-01), where its
  uncertainty is smallest for upcoming semesters; its error, the period error
  and any 3-sigma period aliases go in the extra columns. Circular orbits
  (e = 0).

Usage:
	python targetlist_maker_tess_blind.py
	python targetlist_maker_tess_blind.py --tiers 1 --name blind_UHJ_tier1
	python targetlist_maker_tess_blind.py --ranked ~/Research/TESS/output/blind/candidates_ranked.csv
"""
import argparse
import csv
import os
import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs

warnings.filterwarnings('ignore')

DEFAULT_RANKED = os.path.expanduser('~/Research/TESS/output/blind/candidates_ranked_filtered.csv')
G, MSUN, AU = 6.67430e-11, 1.98847e30, 1.495978707e11
RSUN, RJUP = 6.957e8, 7.1492e7
H_PLANCK, C_LIGHT, K_B = 6.62607015e-34, 2.99792458e8, 1.380649e-23
BANDS = {'H': (1.49e-6, 1.80e-6), 'K': (1.96e-6, 2.48e-6)}  # IGRINS
BOND_ALBEDO, REDISTRIBUTION = 0.1, 0.2


def band_planck(teff, band):
	"""Photon-weighted blackbody flux averaged over a top-hat band (arbitrary units)."""
	wl = np.linspace(*BANDS[band], 200)
	x = H_PLANCK * C_LIGHT / (wl * K_B * teff)
	return np.trapz(wl / (wl ** 5 * np.expm1(np.minimum(x, 700))), wl)


def dayside_contrast(teff, rstar, a_rs, rp_rjup, band):
	"""Planet/star dayside flux ratio in `band`, Cowan & Agol (2011) T_day."""
	t0 = teff / np.sqrt(a_rs)
	tday = t0 * (1 - BOND_ALBEDO) ** 0.25 * (2 / 3. - 5 * REDISTRIBUTION / 12.) ** 0.25
	rp_rs = rp_rjup * RJUP / (rstar * RSUN)
	return rp_rs ** 2 * band_planck(tday, band) / band_planck(teff, band), tday

HEADER = ['Name', 'RA', 'Dec', 'Period (day)', 'a (AU)', 'T transit/inf conj (JD)', 'Obs status',
		  'Eccentric?', 'Kp max (km/s)', 'Transiting?', 'Kmag (<9)', 'e', 'omega', 'inc',
		  'TIC', 'Tier', 'TESS semi-amp (ppm)', 'Implied Rp (RJ)', 'T err (d)', 'P err (d)',
		  'Period aliases (d)', 'Sp type', 'Vmag', 'Teff', 'Notes',
		  'Rstar (Rsun)', 'a/Rstar', 'Tday (K)', 'H contrast (ppm)', 'K contrast (ppm)']


def clean(text):
	"""No commas inside fields (the planners split on commas)."""
	return ' '.join(str(text).replace(',', ';').split())


def query_tic(tic_ids):
	"""One batched TIC query: {TIC ID: row}."""
	result = Catalogs.query_criteria(catalog='Tic', ID=[int(t) for t in tic_ids])
	return {int(row['ID']): row for row in result}


def main():
	p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	p.add_argument('--ranked', default=DEFAULT_RANKED, help='Ranked candidate CSV from tess_blind_rank.py')
	p.add_argument('--tiers', nargs='+', default=['1', '2'], help='Tiers to include (default: 1 2)')
	p.add_argument('--name', default='blind_UHJ_candidates', help='Output: targetlists/<name>_targetlist.csv')
	p.add_argument('--status', default='No data', help='Obs status column value')
	p.add_argument('--rp', type=float, default=1.5, help='Planet radius for the expected contrast [R_J]')
	args = p.parse_args()

	rows = [r for r in csv.DictReader(open(args.ranked)) if r['tier'] in args.tiers]
	seen, keep = set(), []
	for r in rows:  # one entry per star (the best-ranked signal)
		if r['tic'] not in seen:
			seen.add(r['tic'])
			keep.append(r)
	print(f'{len(keep)} stars in tiers {" ".join(args.tiers)} from {args.ranked}')
	tic = query_tic([r['tic'] for r in keep])

	out_rows = []
	for r in keep:
		t = tic.get(int(r['tic']))
		if t is None:
			print('Warning: TIC', r['tic'], 'not found; skipped')
			continue
		coord = SkyCoord(ra=float(t['ra']) * u.deg, dec=float(t['dec']) * u.deg)
		ra = coord.ra.to_string(unit=u.hourangle, sep=':', precision=2, pad=True)
		dec = coord.dec.to_string(unit=u.deg, sep=':', precision=1, alwayssign=True, pad=True)
		mass = float(r['M']) if r['M'] else float(t['mass'])
		per = float(r['P_d'])
		a = (G * mass * MSUN * (per * 86400.) ** 2 / (4 * np.pi ** 2)) ** (1 / 3.)
		kp = 2 * np.pi * a / (per * 86400.) / 1e3
		name = ' '.join(r['name'].split()) or f'TIC {r["tic"]}'
		kmag = float(t['Kmag']) if np.isfinite(float(t['Kmag'])) else ''
		rstar, teff = float(r['R']), float(r['teff'])
		a_rs = a / (rstar * RSUN)
		c_h, tday = dayside_contrast(teff, rstar, a_rs, args.rp, 'H')
		c_k, _ = dayside_contrast(teff, rstar, a_rs, args.rp, 'K')
		notes = r['why'] if r['why'] != 'quiet star, coherent, UHJ-like amplitude' else ''
		if r['flags']:
			notes = (notes + '; ' if notes else '') + r['flags']
		out_rows.append([
			clean(name), ra, dec, f'{per:.6f}', f'{a / AU:.5f}', f'{float(r["T_conj_followup_bjd_tdb"]):.4f}',
			args.status, 'No', f'{kp:.1f}', 'N', f'{kmag:.2f}' if kmag != '' else '', 0, 0, 0,
			r['tic'], r['tier'], f'{float(r["A_ppm"]):.0f}', f'{float(r["implied_rp_rjup"]):.2f}',
			f'{float(r["T_conj_followup_err_d"]):.3f}', f'{float(r["P_err_d"]):.6f}',
			clean(r['aliases_d']), clean(r['sp_type'] or ''), f'{float(r["V"]):.2f}',
			f'{float(r["teff"]):.0f}', clean(notes),
			f'{rstar:.3f}', f'{a_rs:.3f}', f'{tday:.0f}', f'{c_h * 1e6:.0f}', f'{c_k * 1e6:.0f}'])

	os.makedirs('targetlists', exist_ok=True)
	path = os.path.join('targetlists', args.name + '_targetlist.csv')
	with open(path, 'w', newline='') as fh:
		w = csv.writer(fh)
		w.writerow(HEADER)
		w.writerows(out_rows)
	print(f'Saved {len(out_rows)} targets to {path}')


if __name__ == '__main__':
	main()
