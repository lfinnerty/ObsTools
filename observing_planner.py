import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_body
import csv
import sys


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
	site = sys.argv[1]
	if site.upper() in ['KECK', 'GEMINI', 'MKO']:
		site = 'Keck'
		sitestr = 'MKO'
	elif site.upper() in ['MAGELLAN', 'LCO']:
		site = 'Las Campanas Observatory'
		sitestr = 'LCO'
	inst = sys.argv[2]

	### Input: RA, Dec, orbit parameters (a, T, t0), date
	### Output: Visibility, delta Kp

	### load objects and dates here
	### A semester
	# hstdates = ['2026-04-05']
	hstdates =  ['2026-02-01', '2026-02-02', '2026-02-03','2026-02-04','2026-02-05',
				'2026-02-06','2026-02-07','2026-02-08','2026-02-09','2026-02-10',
				'2026-02-11', '2026-02-12', '2026-02-13','2026-02-14','2026-02-15',
				'2026-02-16','2026-02-17','2026-02-18','2026-02-19','2026-02-20',
				'2026-02-21', '2026-02-22', '2026-02-23','2026-02-24','2026-02-25',
				'2026-02-26','2026-02-27','2026-02-28',
				'2026-03-01', '2026-03-02', '2026-03-03','2026-03-04','2026-03-05',
				'2026-03-06','2026-03-07','2026-03-08','2026-03-09','2026-03-10',
				'2026-03-11', '2026-03-12', '2026-03-13','2026-03-14','2026-03-15',
				'2026-03-16','2026-03-17','2026-03-18','2026-03-19','2026-03-20',
				'2026-03-21', '2026-03-22', '2026-03-23','2026-03-24','2026-03-25',
				'2026-03-26','2026-03-27','2026-03-28','2026-03-29','2026-03-30', '2026-03-31',
				'2026-04-01', '2026-04-02', '2026-04-03','2026-04-04','2026-04-04', '2026-04-05',
				'2026-04-06','2026-04-07','2026-04-08','2026-04-09','2026-04-10',
				'2026-04-11', '2026-04-12', '2026-04-13','2026-04-14','2026-04-15',
				'2026-04-16','2026-04-17','2026-04-18','2026-04-19','2026-04-20',
				'2026-04-21', '2026-04-22', '2026-04-23','2026-04-24','2026-04-25',
				'2026-04-26','2026-04-27','2026-04-28','2026-04-29','2026-04-30',
				'2026-05-01', '2026-05-02', '2026-05-03','2026-05-04','2026-05-05',
				'2026-05-06','2026-05-07','2026-05-08','2026-05-09','2026-05-10',
				'2026-05-11', '2026-05-12', '2026-05-13','2026-05-14','2026-05-15',
				'2026-05-16','2026-05-17','2026-05-18','2026-05-19','2026-05-20',
				'2026-05-21', '2026-05-22', '2026-05-23','2026-05-24','2026-05-25',
				'2026-05-26','2026-05-27','2026-05-28','2026-05-29','2026-05-30', '2026-05-31',
				'2026-06-01', '2026-06-02', '2026-06-03','2026-06-04','2026-06-05',
				'2026-06-06','2026-06-07','2026-06-08','2026-06-09','2026-06-10',
				'2026-06-11', '2026-06-12', '2026-06-13','2026-06-14','2026-06-15',
				'2026-06-16','2026-06-17','2026-06-18','2026-06-19','2026-06-20',
				'2026-06-21', '2026-06-22', '2026-06-23','2026-06-24','2026-06-25',
				'2026-06-26','2026-06-27','2026-06-28','2026-06-29','2026-06-30',
				'2026-07-01', '2026-07-02', '2026-07-03','2026-07-04','2026-07-05',
				'2026-07-06','2026-07-07','2026-07-08','2026-07-09','2026-07-10',
				'2026-07-11', '2026-07-12', '2026-07-13','2026-07-14','2026-07-15',
				'2026-07-16','2026-07-17','2026-07-18','2026-07-19','2026-07-20',
				'2026-07-21', '2026-07-22', '2026-07-23','2026-07-24','2026-07-25',
				'2026-07-26','2026-07-27','2026-07-28','2026-07-29','2026-07-30','2026-07-31']
	### B semester
	# hstdates = ['2026-01-21', '2026-01-22', '2026-01-23','2026-01-24','2026-01-25',
	# 			'2026-01-26','2026-01-27','2026-01-28','2026-01-29','2026-01-30','2026-01-31']

	# ['2025-12-01', '2025-12-02', '2025-12-03','2025-12-04','2025-12-05',
	# 			'2025-12-06','2025-12-07','2025-12-08','2025-12-09','2025-12-10',
	# 			'2025-12-11', '2025-12-12', '2025-12-13','2025-12-14','2025-12-15',
	# 			'2025-12-16','2025-12-17','2025-12-18','2025-12-19','2025-12-20',
	# 			'2025-12-21', '2025-12-22', '2025-12-23','2025-12-24','2025-12-25',
	# 			'2025-12-26','2025-12-27','2025-12-28','2025-12-29','2025-12-30','2025-12-31',
	# 			'2026-01-01', '2026-01-02', '2026-01-03','2026-01-04','2026-01-05',
	# 			'2026-01-06','2026-01-07','2026-01-08','2026-01-09','2026-01-10',
	# 			'2026-01-11', '2026-01-12', '2026-01-13','2026-01-14','2026-01-15',
	# 			'2026-01-16','2026-01-17','2026-01-18','2026-01-19','2026-01-20',
	# 			]
	# hstdates = ['2025-08-01', '2025-08-02', '2025-08-03','2025-08-04','2025-08-05',
	# 			'2025-08-06','2025-08-07','2025-08-08','2025-08-09','2025-08-10',
	# 			'2025-08-11', '2025-08-12', '2025-08-13','2025-08-14','2025-08-15',
	# 			'2025-08-16','2025-08-17','2025-08-18','2025-08-19','2025-08-20',
	# 			'2025-08-21', '2025-08-22', '2025-08-23','2025-08-24','2025-08-25',
	# 			'2025-08-26','2025-08-27','2025-08-28','2025-08-29','2025-08-30','2025-08-31',
	# 			'2025-09-01', '2025-09-02', '2025-09-03','2025-09-04','2025-09-05',
	# 			'2025-09-06','2025-09-07','2025-09-08','2025-09-09','2025-09-10',
	# 			'2025-09-11', '2025-09-12', '2025-09-13','2025-09-14','2025-09-15',
	# 			'2025-09-16','2025-09-17','2025-09-18','2025-09-19','2025-09-20',
	# 			'2025-09-21', '2025-09-22', '2025-09-23','2025-09-24','2025-09-25',
	# 			'2025-09-26','2025-09-27','2025-09-28','2025-09-29','2025-09-30',
	# 			'2025-10-01', '2025-10-02', '2025-10-03','2025-10-04','2025-10-05',
	# 			'2025-10-06','2025-10-07','2025-10-08','2025-10-09','2025-10-10',
	# 			'2025-10-11', '2025-10-12', '2025-10-13','2025-10-14','2025-10-15',
	# 			'2025-10-16','2025-10-17','2025-10-18','2025-10-19','2025-10-20',
	# 			'2025-10-21', '2025-10-22', '2025-10-23','2025-10-24','2025-10-25',
	# 			'2025-10-26','2025-10-27','2025-10-28','2025-10-29','2025-10-30', '2025-10-31',
	# 			'2025-11-01', '2025-11-02', '2025-11-03','2025-11-04','2025-11-05',
	# 			'2025-11-06','2025-11-07','2025-11-08','2025-11-09','2025-11-10',
	# 			'2025-11-11', '2025-11-12', '2025-11-13','2025-11-14','2025-11-15',
	# 			'2025-11-16','2025-11-17','2025-11-18','2025-11-19','2025-11-20',
	# 			'2025-11-21', '2025-11-22', '2025-11-23','2025-11-24','2025-11-25',
	# 			'2025-11-26','2025-11-27','2025-11-28','2025-11-29','2025-11-30',]
	# 			'2025-12-01', '2025-12-02', '2025-12-03','2025-12-04','2025-12-05',
	# 			'2025-12-06','2025-12-07','2025-12-08','2025-12-09','2025-12-10',
	# 			'2025-12-11', '2025-12-12', '2025-12-13','2025-12-14','2025-12-15',
	# 			'2025-12-16','2025-12-17','2025-12-18','2025-12-19','2025-12-20',
	# 			'2025-12-21', '2025-12-22', '2025-12-23','2025-12-24','2025-12-25',
	# 			'2025-12-26','2025-12-27','2025-12-28','2025-12-29','2025-12-30','2025-12-31',
	# 			'2025-01-01', '2025-01-02', '2025-01-03','2025-01-04','2025-01-05',
	# 			'2025-01-06','2025-01-07','2025-01-08','2025-01-09','2025-01-10',
	# 			'2025-01-11', '2025-01-12', '2025-01-13','2025-01-14','2025-01-15',
	# 			'2025-01-16','2025-01-17','2025-01-18','2025-01-19','2025-01-20',
	# 			'2025-01-21', '2025-01-22', '2025-01-23','2025-01-24','2025-01-25',
	# 			'2025-01-26','2025-01-27','2025-01-28','2025-01-29','2025-01-30','2025-01-31']
	# hstdates = ['2025-12-07']#,'2025-05-20','2025-05-21','2025-05-22','2025-05-23','2025-05-24','2025-05-25']

	dates = []
	for date in hstdates:
		dates.append(date+'T23:59:59')
	if site == 'Las Campanas Observatory':
		hours = np.linspace(-6,18,100)
	elif site == 'Keck':
		hours = np.linspace(-3,21,100)
	uttimes = hours*u.hour
	for i in range(len(dates)):
		# print(dates[i])
		dates[i] = Time(dates[i], format='isot', scale='utc')+uttimes
	

	objects = np.genfromtxt(inst+'_targetlist.csv', delimiter=',', skip_header=1, dtype=None, encoding=None)


	site = EarthLocation.of_site(site)
	coords = []
	for j in range(len(dates)):
		fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16,12))
		altazframe = AltAz(obstime=dates[j], location=site)
		sunaltaz = get_sun(dates[j]).transform_to(altazframe)
		moonaltaz = get_body('moon',dates[j]).transform_to(altazframe)
		for i in range(len(objects)):
			coord = SkyCoord(ra=objects[i][1], dec=objects[i][2], unit=(u.hourangle, u.deg))
			# print(coord)
			period = objects[i][3]
			kp = objects[i][8] 
			t0 = objects[i][5]

			if objects[i][7] == 'No':
				phase0 = (dates[j].jd - t0)/period
				phase0-=int(phase0[0])
				phis = np.sin(phase0*2*np.pi)
				kps = kp*np.sin(phis)
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
					delta_Kp_min = 50.
				else:
					delta_Kp_min = 40.
					### Set eclipse to 0 for duration
					# transit_kp = kp*np.sin(Rs/a)
					# kps[np.abs(kps)<transit_kp] = 0.

				kpclean = kps[~np.isnan(kps)]
				### Check Kp decreases (i.e. between 0.25 and 0.75 in phase)
				# kp_decreasing = kpclean[-1] < kpclean[0]
				kp_decreasing = True


				if (np.abs(kpmax - kpmin) > delta_Kp_min and kp_decreasing):
					ax[0].plot(uttimes.value[objaltaz.alt>0*u.deg], objaltaz.secz[objaltaz.alt>0*u.deg], label=objects[i][0])
					### Figure out when it's actually above site limits
					abovehorizon = objaltaz.alt.value>25
					rising = np.diff(objaltaz.alt.value, prepend=0) > 0
					overlimit = objaltaz.alt.value>38.
					inds = np.isfinite(kps) & abovehorizon & (rising | overlimit) 
					ax[1].plot(uttimes.value[inds], kps[inds], label=objects[i][0])
					ax[1].text(uttimes.value[inds][-1], kps[inds][-1], objects[i][0])
					# else:
					# 	ax[1].text(uttimes.value[inds][-1], kpclean[-1], objects[i][0])
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
		ax[0].axhline(1./np.cos((90-18)*np.pi/180.), color='r', linestyle='--', label='Rising limit')
		ax[0].axhline(1./np.cos((90-38)*np.pi/180.), color='r', linestyle='-.', label='Setting limit')
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
		ax[1].set_ylabel(r'$v_{pl}$ [km/s]')
		ax[0].set_xlabel('UT time [hr]')
		ax[1].set_xlabel('UT time [hr]')
		ax[0].set_title('HST date: '+hstdates[j])
		plt.savefig('plots/'+hstdates[j]+'_'+sitestr+'.png', bbox_inches='tight')
		plt.close()
		# plt.show()


