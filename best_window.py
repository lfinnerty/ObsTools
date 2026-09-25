import numpy as np
import argparse
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

	objects = np.atleast_1d(np.genfromtxt('targetlists/'+args.instrument+'_targetlist.csv', delimiter=',', skip_header=1, dtype=None, encoding=None))

	rows = []
	for hstdate in generate_hstdates(args.start_date, args.end_date):
		times = Time(hstdate+'T23:59:59', format='isot', scale='utc')+uttimes
		altazframe = AltAz(obstime=times, location=location)
		sunalt = get_sun(times).transform_to(altazframe).alt.deg
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
			period = obj[3]
			t0 = obj[5]
			kp = obj[8]
			ecl_halfwidth = 0.5*args.duration/24./period

			coord = SkyCoord(ra=obj[1], dec=obj[2], unit=(u.hourangle, u.deg))
			alt = coord.transform_to(altazframe).alt.deg
			phase = ((times.jd - t0)/period) % 1
			vpl = kp*np.sin(2*np.pi*phase)
			### Lambertian phase curve (edge-on orbit), normalized to 1 at full phase
			alpha = np.pi*np.abs(1-2*phase)
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
						'crosses_eclipse':np.any(in_eclipse & night & (alt > args.alt_min))})

	### Relative planet S/N
	for r in rows:
		r['snr'] = np.sqrt(r['eff_hrs']*r['dv'])
	snr_max = max([r['snr'] for r in rows]) if rows else 1.
	for r in rows:
		r['snr'] /= snr_max

	### Bright time first, then by planet S/N
	rows.sort(key=lambda r: (not r['bright'], -r['snr']))
	if args.top is not None:
		rows = rows[:args.top]

	print(f"{'rank':>4} {'HST date':10} {'target':14} {'moon':6} {'illum':>5} {'hrs':>5} {'UT start':>8} {'UT end':>7} {'ph0':>5} {'ph1':>5} {'v0':>6} {'v1':>6} {'dv':>5} {'f_pl':>5} {'S/N':>5} {'ecl':>4}")
	for i, r in enumerate(rows):
		print(f"{i+1:4d} {r['date']:10} {r['name'][:14]:14} {'bright' if r['bright'] else 'dark':6} {r['illum']:5.2f} {r['hrs']:5.2f} "
			f"{r['ut0']:8.2f} {r['ut1']:7.2f} {r['ph0']:5.2f} {r['ph1']:5.2f} {r['v0']:6.0f} {r['v1']:6.0f} {r['dv']:5.0f} {r['fpl']:5.2f} {r['snr']:5.2f} {'Y' if r['crosses_eclipse'] else '':>4}")
