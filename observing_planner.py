import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_body
import csv
import sys
import argparse
from datetime import datetime, timedelta


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


def solve_kepler(mobs, e, order=100):
	Eobs = np.copy(mobs)
	for i in range(1,order+1):
		Eobs = mobs + e*np.sin(Eobs)
	return Eobs

def calc_transit_kps(a, Rs, kp0):
	# dphi = np.sin(Rs/a)
	kpmin = -kp0*np.sin(Rs/a)
	kpmax = kp0*np.sin(Rs/a)

if __name__ == '__main__':
	# print(EarthLocation.get_site_names(refresh_cache=True))
	parser = argparse.ArgumentParser()
	parser.add_argument('site')
	parser.add_argument('instrument')
	parser.add_argument('--date', help='single HST date (YYYY-MM-DD)')
	parser.add_argument('--start-date', help='first HST date (YYYY-MM-DD)')
	parser.add_argument('--end-date', help='last HST date (YYYY-MM-DD)')
	parser.add_argument(
		'--require-conjunction',
		action='store_true',
		help='only plot targets whose observed velocity crosses conjunction',
	)
	args = parser.parse_args()
	if args.date is not None and (args.start_date is not None or args.end_date is not None):
		parser.error('--date cannot be combined with --start-date or --end-date')
	if (args.start_date is None) != (args.end_date is None):
		parser.error('--start-date and --end-date must be provided together')
	if args.date is None and args.start_date is None:
		parser.error('provide --date or both --start-date and --end-date')

	site = args.site
	if site.upper() in ['KECK', 'GEMINI-N', 'MKO']:
		site = 'Keck'
		sitestr = 'MKO'
	elif site.upper() in ['MAGELLAN', 'LCO']:
		site = 'Las Campanas Observatory'
		sitestr = 'LCO'
	elif site.upper() in ['EELT', 'VLT', 'PARANAL']:
		site = 'Paranal'
		sitestr = 'EELT'
	elif site.upper() in ['GEMINI-S', 'SOAR', 'CTIO' ]:
		site = 'gemini_south'
		sitestr = 'Gem-S'
	inst = args.instrument

	### Input: RA, Dec, orbit parameters (a, T, t0), date
	### Output: Visibility, delta Kp

	### hardcode dates if you want them
	# hstdates = ['2027-06-24']
	# hstdates = ['2027-01-21', '2027-01-12']
	if args.date is not None:
		hstdates = generate_hstdates(args.date, args.date)
	else:
		hstdates = generate_hstdates(args.start_date, args.end_date)

	dates = []
	for date in hstdates:
		dates.append(date+'T23:59:59')
	if site  in ['Las Campanas Observatory', 'Paranal', 'gemini_south'] :
		hours = np.linspace(-6,18,100)
	elif site in  ['Keck']:
		hours = np.linspace(-3,21,100)
	uttimes = hours*u.hour
	for i in range(len(dates)):
		# print(dates[i])
		dates[i] = Time(dates[i], format='isot', scale='utc')+uttimes
	

	objects = np.atleast_1d(np.genfromtxt('targetlists/'+inst+'_targetlist.csv', delimiter=',', skip_header=1, dtype=None, encoding=None))


	site = EarthLocation.of_site(site)
	coords = []
	for j in range(len(dates)):
		fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16,12))
		altazframe = AltAz(obstime=dates[j], location=site)
		sunaltaz = get_sun(dates[j]).transform_to(altazframe)
		moonaltaz = get_body('moon',dates[j]).transform_to(altazframe)
		nobj = 0
		for i in range(len(objects)):
			coord = SkyCoord(ra=objects[i][1], dec=objects[i][2], unit=(u.hourangle, u.deg))
			# print(coord)
			period = objects[i][3]
			kp = objects[i][8] 
			t0 = objects[i][5]

			if objects[i][7] in ['No', 'N']:
				phase0 = (dates[j].jd - t0)/period
				# print(phase0)
				phase0-=int(phase0[0])
				# print(phase0)
				kps = kp*np.sin(phase0*2*np.pi)
				# print(kps)
			else:
				omega = objects[i][12]*np.pi/180.
				phase0 = (dates[j].jd - t0)/period
				phase0-=int(phase0[0])
				mobs = 2*np.pi*phase0
				ecc = objects[i][11]
				Eobs = solve_kepler(mobs, ecc, order=100)
				beta = ecc/(1+np.sqrt(1-ecc**2))
				fobs = Eobs+2*np.arctan(beta*np.sin(Eobs)/(1-beta*np.cos(Eobs)))
				kps = -kp*(np.cos(fobs+omega)+ecc*np.cos(omega))

			### Now plot object altaz and Kp versus time
			objaltaz = coord.transform_to(AltAz(obstime=dates[j], location=site))
			### Want it above horizon during the night
			if np.sum(objaltaz.alt[sunaltaz.alt < -0*u.deg] > 30*u.deg) > 4:

				### Restrict to nighttime
				kps[sunaltaz.alt > -0*u.deg] = np.nan
				kpmin = np.nanmin(kps[objaltaz.alt>30*u.deg])
				kpmax = np.nanmax(kps[objaltaz.alt>30*u.deg])
				
				### Minimum Kp shift - bigger if non-transiting. Also figure out when transit is to mask it
				if objects[i][9] == 'N':
					delta_Kp_min = 30.
				else:
					delta_Kp_min = 30.
					### Set eclipse to 0 for duration
					# transit_kp = kp*np.sin(Rs/a)
					# kps[np.abs(kps)<transit_kp] = 0.

				kpclean = kps[~np.isnan(kps)]
				### Check Kp decreases (i.e. between 0.25 and 0.75 in phase)
				kp_decreasing = kpclean[-1] < kpclean[0]


				if (np.abs(kpmax - kpmin)) > delta_Kp_min and kp_decreasing:
					# if objects[i][0] in ['WASP-121']:
					# print('Making plot')	
					ax[0].plot(uttimes.value[objaltaz.alt>0*u.deg], objaltaz.secz[objaltaz.alt>0*u.deg], label=objects[i][0])
					### Figure out when it's actually above site limits
					abovehorizon = objaltaz.alt.value>15
					rising = np.diff(objaltaz.alt.value, prepend=0) > 0
					overlimit = objaltaz.alt.value>15#38.
					inds = np.isfinite(kps) & abovehorizon & (rising | overlimit) 
					if (not args.require_conjunction):
						ax[1].plot(uttimes.value[inds], kps[inds], label=objects[i][0])
						ax[1].text(uttimes.value[inds][-1], kps[inds][-1], objects[i][0])
						nobj+=1
					elif args.require_conjunction:
						if kps[inds][0] > 0 and kps[inds][-1] < 0:
							ax[1].plot(uttimes.value[inds], kps[inds], label=objects[i][0])
							ax[1].text(uttimes.value[inds][-1], kps[inds][-1], objects[i][0])
							nobj+=1
					if nobj == 1:
						print('Making plot for', hstdates[j], 'at', sitestr)
		### overlay sun
		ax[0].plot(uttimes.value, sunaltaz.secz, color='y', label='Sun')
		### And moon
		# mooon = moonaltaz.secz
		# moon[moon<1] = np.nan
		ax[0].plot(uttimes.value, moonaltaz.secz, color='r', label='Moon')
		ax[0].fill_between(uttimes.value, 0, 90, sunaltaz.alt < -0*u.deg, color='0.8', zorder=0) 
		ax[0].fill_between(uttimes.value, 0, 90, sunaltaz.alt < -18*u.deg, color='0.6', zorder=0) 
		ax[1].fill_between(uttimes.value, -300,300, sunaltaz.alt < -0*u.deg, color='0.8', zorder=0) 
		ax[1].fill_between(uttimes.value, -300,300, sunaltaz.alt < -18*u.deg, color='0.6', zorder=0) 
		ax[1].axhline(0,color='k',linestyle='--')
		# ax[0].axhline(1./np.cos((90-18)*np.pi/180.), color='r', linestyle='--', label='Rising limit')
		# ax[0].axhline(1./np.cos((90-38)*np.pi/180.), color='r', linestyle='-.', label='Setting limit')
		ax[0].xaxis.set_major_locator(ticker.FixedLocator(np.arange(hours[0],hours[-1]+0.01,1)))
		ax[1].xaxis.set_major_locator(ticker.FixedLocator(np.arange(hours[0],hours[-1]+0.01,1)))

		ax[0].grid(visible=True)
		ax[1].grid(visible=True)
		ax[0].legend()
		# ax[1].legend()

		sunset = np.min(uttimes.value[sunaltaz.alt<-0*u.deg])-1
		sunrise = np.max(uttimes.value[sunaltaz.alt<-0*u.deg])+1

		ax[0].set_ylim(0, 6)
		ax[1].set_ylim(-200,200)
		ax[0].set_xlim(sunset, sunrise)
		ax[1].set_xlim(sunset,sunrise)
		ax[0].set_ylabel('Airmass')
		ax[0].set_ylim([1.0,3.0])
		ax[1].set_ylabel(r'$v_{pl}$ [km/s]')
		ax[0].set_xlabel('UT time [hr]')
		ax[1].set_xlabel('UT time [hr]')
		ax[0].set_title('HST date: '+hstdates[j])
		if nobj>0:
			plt.savefig('plots/'+hstdates[j]+'_'+sitestr+'.png', bbox_inches='tight')
		plt.close()
		# plt.show()


