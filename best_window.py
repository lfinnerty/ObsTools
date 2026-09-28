import numpy as np
import argparse
import csv
import warnings
from datetime import datetime, timedelta
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_body

warnings.filterwarnings('ignore')

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
	parser.add_argument('--csv', default=None, help='save the ranked windows to this CSV file')
	parser.add_argument('--weight-kmag', action='store_true', help='weight the S/N by stellar K-band photon flux, 10^(-0.2 (Kmag - 6)) (column 10); targets without Kmag are skipped')
	parser.add_argument('--weight-contrast', action='store_true', help="weight the S/N by the planet's expected dayside contrast (target-list column named by --contrast-column)")
	parser.add_argument('--contrast-column', default='K contrast (ppm)', help='header name of the contrast column used by --weight-contrast')
	parser.add_argument('--snr-min', type=float, default=0.3, help='drop windows with relative S/N below this (0 keeps all)')
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
			usable_hrs = np.sum(usable)*dt_hr
			### Effective hours, weighting each exposure by f^2
			eff_hrs = np.sum(fpl[usable]**2)*dt_hr
			delta_v = np.max(vpl[usable]) - np.min(vpl[usable])
			if delta_v < args.delta_v_min:
				continue

			ut_usable = hours[usable]
			rows.append({'date':hstdate, 'name':obj[0], 'bright':bright, 'illum':illum,
						'hrs':usable_hrs, 'eff_hrs':eff_hrs, 'fpl':np.mean(fpl[usable]), 'ut0':ut_usable[0], 'ut1':ut_usable[-1],
						'ph0':phase[usable][0], 'ph1':phase[usable][-1],
						'v0':vpl[usable][0], 'v1':vpl[usable][-1], 'dv':delta_v,
						'crosses_eclipse':np.any(in_eclipse & night & (alt > args.alt_min)), 'weight':weight})

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

	print(f"{'rank':>4} {'HST date':10} {'target':14} {'moon':6} {'illum':>5} {'hrs':>5} {'UT start':>8} {'UT end':>7} {'ph0':>5} {'ph1':>5} {'v0':>6} {'v1':>6} {'dv':>5} {'f_pl':>5} {'S/N':>5} {'ecl':>4}")
	for i, r in enumerate(rows):
		print(f"{i+1:4d} {r['date']:10} {r['name'][:14]:14} {'bright' if r['bright'] else 'dark':6} {r['illum']:5.2f} {r['hrs']:5.2f} "
			f"{r['ut0']:8.2f} {r['ut1']:7.2f} {r['ph0']:5.2f} {r['ph1']:5.2f} {r['v0']:6.0f} {r['v1']:6.0f} {r['dv']:5.0f} {r['fpl']:5.2f} {r['snr']:5.2f} {'Y' if r['crosses_eclipse'] else '':>4}")

	if args.csv is not None:
		with open(args.csv, 'w', newline='') as f:
			writer = csv.writer(f)
			writer.writerow(['rank', 'hst_date', 'target', 'moon', 'illum', 'usable_hrs', 'ut_start', 'ut_end',
							'phase_start', 'phase_end', 'v_start', 'v_end', 'delta_v', 'f_pl', 'rel_snr', 'crosses_eclipse'])
			for i, r in enumerate(rows):
				writer.writerow([i+1, r['date'], r['name'], 'bright' if r['bright'] else 'dark', f"{r['illum']:.2f}", f"{r['hrs']:.2f}",
								f"{r['ut0']:.2f}", f"{r['ut1']:.2f}", f"{r['ph0']:.3f}", f"{r['ph1']:.3f}", f"{r['v0']:.1f}", f"{r['v1']:.1f}",
								f"{r['dv']:.1f}", f"{r['fpl']:.3f}", f"{r['snr']:.3f}", 'Y' if r['crosses_eclipse'] else 'N'])
		print('Saved', args.csv)
