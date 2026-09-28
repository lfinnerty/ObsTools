import numpy as np
import argparse
import csv
import os
import re
import sys
import warnings
from datetime import datetime, timedelta
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_body
from numpy.lib.stride_tricks import sliding_window_view

warnings.filterwarnings('ignore')

OUTPUT_DIR = 'output/best_windows'

### Rank nights for dayside emission HRCCS of circular-orbit planets.
### Usable time: target above alt_min, sun below -18 deg (no astronomical twilight),
### planet on dayside half of orbit (phase 0.25-0.75), and outside secondary eclipse.
### Bright-time nights rank above dark-time nights (IR is insensitive to sky brightness);
### within each category, nights are ranked by expected planet S/N.
### Planet S/N scales as (stellar S/N) * sqrt(delta v), with stellar S/N ~ sqrt(usable time).
### Optionally (--weight-kmag, --weight-contrast) it is multiplied by 10^(-0.2 Kmag), for
### photon-noise-limited stars, and by the planet's dayside contrast, so that S/N compares
### targets as well as nights. Windows below --snr-min (default 0.3) are dropped.
### The planet signal in each exposure is scaled by a Lambertian phase curve,
### f = (sin(a) + (pi - a) cos(a))/pi, with phase angle a = pi |1 - 2 phase|
### (a = 0 at secondary eclipse), and a matched-filter sum over exposures gives
### S/N ~ sqrt(sum(f^2 dt) * delta v), normalized to the best night.
### With --max-hours, each night's window is limited to that much wall-clock time (for
### targets that do not need the whole night): the contiguous stretch of at most that
### length, within the usable time, that maximizes sum(f^2 dt) * delta v is used.
### With --airmass-k and/or --seeing-exponent, each exposure is also weighted by the airmass X
### (photon-noise limited, so the weight multiplies S/N^2): atmospheric extinction
### 10^(-0.4 k (X-1)) with k in mag/airmass, and slit losses for seeing FWHM ~ X^0.6 wider
### than the slit, (1/X)^p with p = --seeing-exponent. The 4 h cap uses the same weights.
### --csv files without a directory are written to OUTPUT_DIR, together with a .txt copy
### of the printed table; --split-targets also writes one file per target.


def generate_hstdates(start_date, end_date):
	"""Return ISO-formatted dates from start_date through end_date, inclusive."""
	start = datetime.strptime(start_date, '%Y-%m-%d').date()
	end = datetime.strptime(end_date, '%Y-%m-%d').date()
	if start > end:
		raise ValueError('start date must not be after end date')

	hstdates = []
	current = start
	while current <= end:
		hstdates.append(current.isoformat())
		current += timedelta(days=1)
	return hstdates


def best_subwindow(usable, w2, vpl, nmax):
	"""Mask of the stretch of at most nmax consecutive samples that maximizes
	sum(w2 over usable samples) * (planet velocity range over usable samples),
	where w2 is the per-sample S/N^2 weight (phase curve f^2 x airmass weight)."""
	idx = np.flatnonzero(usable)
	i0, i1 = idx[0], idx[-1]+1
	if i1-i0 <= nmax:
		return usable
	u = usable[i0:i1]
	f2 = np.where(u, w2[i0:i1], 0.)
	v = np.where(u, vpl[i0:i1], np.nan)
	eff = sliding_window_view(f2, nmax).sum(axis=1)
	vw = sliding_window_view(v, nmax)
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', RuntimeWarning)
		dv = np.nan_to_num(np.nanmax(vw, axis=1) - np.nanmin(vw, axis=1))
	k = np.argmax(eff*dv)
	out = np.zeros_like(usable)
	out[i0+k:i0+k+nmax] = u[k:k+nmax]
	return out


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('site')
	parser.add_argument('instrument')
	parser.add_argument('--start-date', required=True, help='first HST date (YYYY-MM-DD)')
	parser.add_argument('--end-date', required=True, help='last HST date (YYYY-MM-DD)')
	parser.add_argument('--duration', type=float, required=True, help='transit/eclipse duration (hours), applied to all targets')
	parser.add_argument('--alt-min', type=float, default=30., help='minimum target altitude (deg)')
	parser.add_argument('--delta-v-min', type=float, default=30., help='minimum planet velocity change over usable time (km/s)')
	parser.add_argument('--bright-threshold', type=float, default=0.5, help='moon illumination at local midnight above which a night is bright time')
	parser.add_argument('--top', type=int, default=None, help='only print the top N nights')
	parser.add_argument('--moon-sep-min', type=float, default=10., help='minimum target-moon separation (deg) while the target is up at night')
	parser.add_argument('--csv', default=None, help=f'save the ranked windows to this CSV file (in {OUTPUT_DIR}/ unless a directory is given) and the printed table to the same name with .txt')
	parser.add_argument('--weight-kmag', action='store_true', help='weight the S/N by stellar K-band photon flux, 10^(-0.2 (Kmag - 6)) (column 10); targets without Kmag are skipped')
	parser.add_argument('--weight-contrast', action='store_true', help="weight the S/N by the planet's expected dayside contrast (target-list column named by --contrast-column)")
	parser.add_argument('--contrast-column', default='K contrast (ppm)', help='header name of the contrast column used by --weight-contrast')
	parser.add_argument('--split-targets', action='store_true', help='with --csv, also write one CSV/.txt per target (<csv name>_<target>.csv), ranked within the target')
	parser.add_argument('--snr-min', type=float, default=0.3, help='drop windows with relative S/N below this (0 keeps all)')
	parser.add_argument('--max-hours', type=float, default=None, help='maximum wall-clock length (hours) of a window: the best stretch of this length within the usable time is used, for targets that do not need the whole night')
	parser.add_argument('--airmass-k', type=float, default=0., help='extinction (mag/airmass) for the airmass weight 10^(-0.4 k (X-1)) on S/N^2 (0: none; ~0.05 in H/K at Maunakea)')
	parser.add_argument('--seeing-exponent', type=float, default=0., help='slit-loss airmass weight X^(-p) on S/N^2, for seeing FWHM ~ X^0.6 wider than the slit (0: none; 0.6 in that limit)')
	parser.add_argument('--inclination', type=float, default=90., help='orbital inclination (deg) applied to all targets: Kp = Kp_max sin(i), and the Lambertian phase angle uses cos(a) = -sin(i) cos(2 pi phase)')
	args = parser.parse_args()

	site = args.site
	if site.upper() in ['KECK', 'GEMINI-N', 'MKO']:
		site = 'Keck'
		hours = np.linspace(-3,21,1441)
	elif site.upper() in ['MAGELLAN', 'LCO']:
		site = 'Las Campanas Observatory'
		hours = np.linspace(-6,18,1441)
	elif site.upper() in ['EELT', 'VLT', 'PARANAL']:
		site = 'Paranal'
		hours = np.linspace(-6,18,1441)
	elif site.upper() in ['GEMINI-S', 'SOAR', 'CTIO' ]:
		site = 'gemini_south'
		hours = np.linspace(-6,18,1441)
	location = EarthLocation.of_site(site)
	dt_hr = hours[1]-hours[0]
	uttimes = hours*u.hour

	tlfile = 'targetlists/'+args.instrument+'_targetlist.csv'
	objects = np.atleast_1d(np.genfromtxt(tlfile, delimiter=',', skip_header=1, dtype=None, encoding=None))
	header = open(tlfile).readline().strip().split(',')
	icontrast = None
	if args.weight_contrast:
		if args.contrast_column not in header:
			raise SystemExit(f"--weight-contrast: no '{args.contrast_column}' column in {tlfile}")
		icontrast = header.index(args.contrast_column)

	### Per-target S/N weights: stellar photon flux (photon-noise limited) and planet contrast
	def target_weight(obj):
		w = 1.
		if args.weight_kmag:
			try:
				kmag = float(obj[10])
			except (TypeError, ValueError):
				kmag = np.nan
			if not np.isfinite(kmag):
				return np.nan
			w *= 10**(-0.2*(kmag-6.))
		if icontrast is not None:
			try:
				w *= float(obj[icontrast])
			except (TypeError, ValueError):
				return np.nan
		return w
	skipped = [obj[0] for obj in objects if not np.isfinite(target_weight(obj))]
	if skipped:
		print('Skipping (missing Kmag or contrast):', ', '.join(skipped))

	rows = []
	for hstdate in generate_hstdates(args.start_date, args.end_date):
		times = Time(hstdate+'T23:59:59', format='isot', scale='utc')+uttimes
		altazframe = AltAz(obstime=times, location=location)
		sunalt = get_sun(times).transform_to(altazframe).alt.deg
		moonaltaz = get_body('moon', times).transform_to(altazframe)
		night = sunalt < -18.

		### Moon illumination at local midnight (middle of the night)
		tmid = times[night][len(times[night])//2] if np.any(night) else times[len(times)//2]
		elong = get_body('moon', tmid).separation(get_sun(tmid)).rad
		illum = (1-np.cos(elong))/2.
		bright = illum >= args.bright_threshold

		for obj in objects:
			if obj[7] not in ['No', 'N']:
				print('Skipping', obj[0], '- eccentric orbits not supported')
				continue
			weight = target_weight(obj)
			if not np.isfinite(weight):
				continue
			period = obj[3]
			t0 = obj[5]
			kp = obj[8]*np.sin(np.radians(args.inclination))
			ecl_halfwidth = 0.5*args.duration/24./period

			coord = SkyCoord(ra=obj[1], dec=obj[2], unit=(u.hourangle, u.deg))
			objaltaz = coord.transform_to(altazframe)
			alt = objaltaz.alt.deg
			airmass = 1./np.sin(np.radians(np.clip(alt, 1., 90.)))
			wam = 10**(-0.4*args.airmass_k*(airmass-1.)) * airmass**(-args.seeing_exponent)
			### Skip if the moon comes within moon_sep_min of the target while it's up at night
			upatnight = (sunalt < 0.) & (alt > 0.)
			if np.any(objaltaz.separation(moonaltaz).deg[upatnight] < args.moon_sep_min):
				continue
			phase = ((times.jd - t0)/period) % 1
			vpl = kp*np.sin(2*np.pi*phase)
			### Lambertian phase curve, normalized to 1 at full phase; phase angle a with
			### cos(a) = -sin(i) cos(2 pi phase) (edge-on: a = pi |1 - 2 phase|)
			alpha = np.arccos(np.clip(-np.sin(np.radians(args.inclination))*np.cos(2*np.pi*phase), -1, 1))
			fpl = (np.sin(alpha) + (np.pi-alpha)*np.cos(alpha))/np.pi

			dayside = (phase > 0.25) & (phase < 0.75)
			in_eclipse = np.abs(phase-0.5) < ecl_halfwidth
			usable = night & (alt > args.alt_min) & dayside & ~in_eclipse
			if not np.any(usable):
				continue
			if args.max_hours is not None:
				usable = best_subwindow(usable, fpl**2*wam, vpl, int(round(args.max_hours/dt_hr))+1)
			usable_hrs = np.sum(usable)*dt_hr
			### Effective hours, weighting each exposure by f^2
			eff_hrs = np.sum(fpl[usable]**2*wam[usable])*dt_hr
			delta_v = np.max(vpl[usable]) - np.min(vpl[usable])
			if delta_v < args.delta_v_min:
				continue

			ut_usable = hours[usable]
			span = (hours >= ut_usable[0]) & (hours <= ut_usable[-1])
			rows.append({'date':hstdate, 'name':obj[0], 'bright':bright, 'illum':illum,
						'hrs':usable_hrs, 'eff_hrs':eff_hrs, 'fpl':np.mean(fpl[usable]), 'ut0':ut_usable[0], 'ut1':ut_usable[-1],
						'X_mean':np.mean(airmass[usable]), 'X_max':np.max(airmass[usable]),
						'ph0':phase[usable][0], 'ph1':phase[usable][-1],
						'v0':vpl[usable][0], 'v1':vpl[usable][-1], 'dv':delta_v,
						'crosses_eclipse':np.any(in_eclipse & night & (alt > args.alt_min) & span), 'weight':weight})

	### Relative planet S/N (x stellar flux and planet contrast weights if requested)
	for r in rows:
		r['snr'] = np.sqrt(r['eff_hrs']*r['dv'])*r['weight']
	snr_max = max([r['snr'] for r in rows]) if rows else 1.
	for r in rows:
		r['snr'] /= snr_max
	n_all = len(rows)
	rows = [r for r in rows if r['snr'] >= args.snr_min]
	print(f'{len(rows)} of {n_all} windows with relative S/N >= {args.snr_min}')

	### Bright time first, then by planet S/N
	rows.sort(key=lambda r: (not r['bright'], -r['snr']))
	if args.top is not None:
		rows = rows[:args.top]

	def table_lines(rows, first):
		lines = [first,
			f"{'rank':>4} {'HST date':10} {'target':14} {'moon':6} {'illum':>5} {'hrs':>5} {'UT start':>8} {'UT end':>7} {'ph0':>5} {'ph1':>5} {'v0':>6} {'v1':>6} {'dv':>5} {'f_pl':>5} {'X':>4} {'Xmax':>4} {'S/N':>5} {'ecl':>4}"]
		for i, r in enumerate(rows):
			lines.append(f"{i+1:4d} {r['date']:10} {r['name'][:14]:14} {'bright' if r['bright'] else 'dark':6} {r['illum']:5.2f} {r['hrs']:5.2f} "
				f"{r['ut0']:8.2f} {r['ut1']:7.2f} {r['ph0']:5.2f} {r['ph1']:5.2f} {r['v0']:6.0f} {r['v1']:6.0f} {r['dv']:5.0f} {r['fpl']:5.2f} {r['X_mean']:4.2f} {r['X_max']:4.2f} {r['snr']:5.2f} {'Y' if r['crosses_eclipse'] else '':>4}")
		return lines

	def write_outputs(path, rows, first):
		"""Ranked windows to path (CSV) and the printed table, with the command, to .txt."""
		with open(os.path.splitext(path)[0]+'.txt', 'w') as f:
			f.write(' '.join(sys.argv)+'\n'+'\n'.join(table_lines(rows, first))+'\n')
		with open(path, 'w', newline='') as f:
			writer = csv.writer(f)
			writer.writerow(['rank', 'hst_date', 'target', 'moon', 'illum', 'usable_hrs', 'ut_start', 'ut_end',
							'phase_start', 'phase_end', 'v_start', 'v_end', 'delta_v', 'f_pl', 'rel_snr', 'crosses_eclipse', 'airmass_mean', 'airmass_max'])
			for i, r in enumerate(rows):
				writer.writerow([i+1, r['date'], r['name'], 'bright' if r['bright'] else 'dark', f"{r['illum']:.2f}", f"{r['hrs']:.2f}",
								f"{r['ut0']:.2f}", f"{r['ut1']:.2f}", f"{r['ph0']:.3f}", f"{r['ph1']:.3f}", f"{r['v0']:.1f}", f"{r['v1']:.1f}",
								f"{r['dv']:.1f}", f"{r['fpl']:.3f}", f"{r['snr']:.3f}", 'Y' if r['crosses_eclipse'] else 'N', f"{r['X_mean']:.3f}", f"{r['X_max']:.3f}"])
		print('Saved', path, 'and', os.path.splitext(path)[0]+'.txt')

	first = f"{len(rows)} of {n_all} windows with relative S/N >= {args.snr_min}"
	print('\n'.join(table_lines(rows, first)[1:]))

	if args.csv is not None:
		if not os.path.dirname(args.csv):
			args.csv = os.path.join(OUTPUT_DIR, args.csv)
		os.makedirs(os.path.dirname(args.csv), exist_ok=True)
		write_outputs(args.csv, rows, first)
		if args.split_targets:
			stem, ext = os.path.splitext(args.csv)
			for name in dict.fromkeys(r['name'] for r in rows):
				sub = [r for r in rows if r['name'] == name]
				write_outputs(f"{stem}_{re.sub(r'[^A-Za-z0-9]', '', name)}{ext}", sub,
							f"{name}: {len(sub)} windows (relative S/N normalized over all targets in the list)")
