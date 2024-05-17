import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import csv


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
	### Input: RA, Dec, orbit parameters (a, T, t0), date
	### Output: Visibility, delta Kp

	### load objects and dates here
	# ### A semester
	hstdates = ['2024-04-01', '2024-04-02', '2024-04-03','2024-04-04','2024-04-04',
				'2024-04-06','2024-04-07','2024-04-08','2024-04-09','2024-04-10',
				'2024-04-11', '2024-04-12', '2024-04-13','2024-04-14','2024-04-15',
				'2024-04-16','2024-04-17','2024-04-18','2024-04-19','2024-04-20',
				'2024-04-21', '2024-04-22', '2024-04-23','2024-04-24','2024-04-25',
				'2024-04-26','2024-04-27','2024-04-28','2024-04-29','2024-04-30']

	# '2024-05-01', '2024-05-02', '2024-05-03','2024-05-04','2024-05-05',
	# 			'2024-05-06','2024-05-07','2024-05-08','2024-05-09','2024-05-10',
	# 			'2024-05-11', '2024-05-12', '2024-05-13','2024-05-14','2024-05-15',
	# 			'2024-05-16','2024-05-17','2024-05-18','2024-05-19','2024-05-20',
	# 			'2024-05-21', '2024-05-22', '2024-05-23','2024-05-24','2024-05-25',
	# 			'2024-05-26','2024-05-27','2024-05-28','2024-05-29','2024-05-30', '2024-05-31',
	# 			'2024-06-01', '2024-06-02', '2024-06-03','2024-06-04','2024-06-05',
	# 			'2024-06-06','2024-06-07','2024-06-08','2024-06-09','2024-06-10',
	# 			'2024-06-11', '2024-06-12', '2024-06-13','2024-06-14','2024-06-15',
	# 			'2024-06-16','2024-06-17','2024-06-18','2024-06-19','2024-06-20',
	# 			'2024-06-21', '2024-06-22', '2024-06-23','2024-06-24','2024-06-25',
	# 			'2024-06-26','2024-06-27','2024-06-28','2024-06-29','2024-06-30',
	# 			'2024-07-01', '2024-07-02', '2024-07-03','2024-07-04','2024-07-05',
	# 			'2024-07-06','2024-07-07','2024-07-08','2024-07-09','2024-07-10',
	# 			'2024-07-11', '2024-07-12', '2024-07-13','2024-07-14','2024-07-15',
	# 			'2024-07-16','2024-07-17','2024-07-18','2024-07-19','2024-07-20',
	# 			'2024-07-21', '2024-07-22', '2024-07-23','2024-07-24','2024-07-25',
	# 			'2024-07-26','2024-07-27','2024-07-28','2024-07-29','2024-07-30','2024-07-31']
	### B semester
	# hstdates = ['2024-08-01', '2024-08-02', '2024-08-03','2024-08-04','2024-08-05',
	# 			'2024-08-06','2024-08-07','2024-08-08','2024-08-09','2024-08-10',
	# 			'2024-08-11', '2024-08-12', '2024-08-13','2024-08-14','2024-08-15',
	# 			'2024-08-16','2024-08-17','2024-08-18','2024-08-19','2024-08-20',
	# 			'2024-08-21', '2024-08-22', '2024-08-23','2024-08-24','2024-08-25',
	# 			'2024-08-26','2024-08-27','2024-08-28','2024-08-29','2024-08-30','2024-08-31',
	# 			'2024-09-01', '2024-09-02', '2024-09-03','2024-09-04','2024-09-05',
	# 			'2024-09-06','2024-09-07','2024-09-08','2024-09-09','2024-09-10',
	# 			'2024-09-11', '2024-09-12', '2024-09-13','2024-09-14','2024-09-15',
	# 			'2024-09-16','2024-09-17','2024-09-18','2024-09-19','2024-09-20',
	# 			'2024-09-21', '2024-09-22', '2024-09-23','2024-09-24','2024-09-25',
	# 			'2024-09-26','2024-09-27','2024-09-28','2024-09-29','2024-09-30',
	# 			'2024-10-01', '2024-10-02', '2024-10-03','2024-10-04','2024-10-05',
	# 			'2024-10-06','2024-10-07','2024-10-08','2024-10-09','2024-10-10',
	# 			'2024-10-11', '2024-10-12', '2024-10-13','2024-10-14','2024-10-15',
	# 			'2024-10-16','2024-10-17','2024-10-18','2024-10-19','2024-10-20',
	# 			'2024-10-21', '2024-10-22', '2024-10-23','2024-10-24','2024-10-25',
	# 			'2024-10-26','2024-10-27','2024-10-28','2024-10-29','2024-10-30', '2024-10-31',
	# 			'2024-11-01', '2024-11-02', '2024-11-03','2024-11-04','2024-11-05',
	# 			'2024-11-06','2024-11-07','2024-11-08','2024-11-09','2024-11-10',
	# 			'2024-11-11', '2024-11-12', '2024-11-13','2024-11-14','2024-11-15',
	# 			'2024-11-16','2024-11-17','2024-11-18','2024-11-19','2024-11-20',
	# 			'2024-11-21', '2024-11-22', '2024-11-23','2024-11-24','2024-11-25',
	# 			'2024-11-26','2024-11-27','2024-11-28','2024-11-29','2024-11-30',
	# 			'2024-12-01', '2024-12-02', '2024-12-03','2024-12-04','2024-12-05',
	# 			'2024-12-06','2024-12-07','2024-12-08','2024-12-09','2024-12-10',
	# 			'2024-12-11', '2024-12-12', '2024-12-13','2024-12-14','2024-12-15',
	# 			'2024-12-16','2024-12-17','2024-12-18','2024-12-19','2024-12-20',
	# 			'2024-12-21', '2024-12-22', '2024-12-23','2024-12-24','2024-12-25',
	# 			'2024-12-26','2024-12-27','2024-12-28','2024-12-29','2024-12-30','2024-12-31',
	# 			'2024-12-01', '2024-12-02', '2024-12-03','2024-12-04','2024-12-05',
	# 			'2025-01-06','2025-01-07','2025-01-08','2025-01-09','2025-01-10',
	# 			'2025-01-11', '2025-01-12', '2025-01-13','2025-01-14','2025-01-15',
	# 			'2025-01-16','2025-01-17','2025-01-18','2025-01-19','2025-01-20',
	# 			'2025-01-21', '2025-01-22', '2025-01-23','2025-01-24','2025-01-25',
	# 			'2025-01-26','2025-01-27','2025-01-28','2025-01-29','2025-01-30','2025-01-31']
	# hstdates = ['2024-05-19','2024-05-20','2024-05-21','2024-05-22','2024-05-23','2024-05-24','2024-05-25']

	dates = []
	for date in hstdates:
		dates.append(date+'T23:59:59')
	uttimes = +np.linspace(3,18,100)*u.hour
	for i in range(len(dates)):
		# print(dates[i])
		dates[i] = Time(dates[i], format='isot', scale='utc')+uttimes
	
	# with open('obsplanning_targetlist.csv') as infile:
	# 	objects = csv.reader(infile, delimiter=',')
	objects = np.genfromtxt('obsplanning_targetlist_2024B.csv', delimiter=',', skip_header=1, dtype=None, encoding=None)
	# print(objects)

	# objects = [['tau Boo', '13:47:17.0', '+17:27:22', 3.31249, 0.046, 2450529.2]]
			   

	keck = EarthLocation.of_site('Keck')
	coords = []
	for j in range(len(dates)):
		fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16,12))
		altazframe = AltAz(obstime=dates[j], location=keck)
		sunaltaz = get_sun(dates[j]).transform_to(altazframe)
		for i in range(len(objects)):
			### Skip eccentric velocities for now
			# if objects[i][7] == 'Yes':
			# 	continue
			# if objects[i][6] in ['KPIC detected', 'KPIC published', 'KPIC tentative', 'KPIC not analyzed']:
			# 	continue
			# if objects[i][0] in ['KELT-7', 'HD93963 b',  'HD 93963 c', 'HD 143105', 'HAT-P-22', 'HIP 86221', 'HD 103720',
			# 					'WASP-189', 'HD 149143', 'GJ 1214', 'MASCARA-1', 'KELT20/M-2', 'HD 149026', 'TOI 1431',
			# 					'TOI 1518', 'HD 209458', 'HD 189733', '51 Peg', 'WASP-33', 'KELT-9']:
			# 					continue
			# if objects[i][0] not in ['HD 189733', '51 Peg', 'tau Boo']:
				# continue
			# print(str(objects[i][1]))
			coord = SkyCoord(ra=objects[i][1], dec=objects[i][2], unit=(u.hourangle, u.deg))
			# print(coord)
			period = objects[i][3]
			kp = objects[i][8] #2*np.pi*objects[i][4]*1.5e8/(period*86400)
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
			objaltaz = coord.transform_to(AltAz(obstime=dates[j], location=keck))
			### Want it above horizon during the night
			if np.sum(objaltaz.alt[sunaltaz.alt < -0*u.deg] > 30*u.deg) > 4:
				### Restrict to nighttime
				kps[sunaltaz.alt > -0*u.deg] = np.nan
				kpmin = np.nanmin(kps[objaltaz.alt>30*u.deg])
				kpmax = np.nanmax(kps[objaltaz.alt>30*u.deg])
				
				### Minimum Kp shift - bigger if non-transiting. Also figure out when transit is to mask it
				if objects[i][9] == 'N':
					delta_Kp_min = 60.
				else:
					delta_Kp_min = 40.
					### Set eclipse to 0 for duration
					# transit_kp = kp*np.sin(Rs/a)
					# kps[np.abs(kps)<transit_kp] = 0.

				kpclean = kps[~np.isnan(kps)]
				### Check Kp decreases (i.e. between 0.25 and 0.75 in phase)
				kp_decreasing = kpclean[-1] < kpclean[0]


				if (np.abs(kpmax - kpmin) > delta_Kp_min and kp_decreasing):# or objects[i][0] in ['WASP-14','HD202772']:
					ax[0].plot(uttimes.value[objaltaz.alt>0*u.deg], objaltaz.secz[objaltaz.alt>0*u.deg], label=objects[i][0])
					### Figure out when it's actually above Keck limits
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
		ax[0].fill_between(uttimes.value, 0, 90, sunaltaz.alt < -0*u.deg, color='0.8', zorder=0) 
		ax[0].fill_between(uttimes.value, 0, 90, sunaltaz.alt < -18*u.deg, color='0.6', zorder=0) 
		ax[1].fill_between(uttimes.value, -300,300, sunaltaz.alt < -0*u.deg, color='0.8', zorder=0) 
		ax[1].fill_between(uttimes.value, -300,300, sunaltaz.alt < -18*u.deg, color='0.6', zorder=0) 
		ax[0].axhline(1./np.cos((90-18)*np.pi/180.), color='r', linestyle='--', label='Rising limit')
		ax[0].axhline(1./np.cos((90-38)*np.pi/180.), color='r', linestyle='-.', label='Setting limit')
		ax[0].xaxis.set_major_locator(ticker.FixedLocator(np.arange(3,18.5,1)))
		ax[1].xaxis.set_major_locator(ticker.FixedLocator(np.arange(3,18.5,1)))

		ax[0].grid(visible=True)
		ax[1].grid(visible=True)
		ax[0].legend()
		# ax[1].legend()
		ax[0].set_ylim(0, 6)
		ax[1].set_ylim(-200,200)
		ax[0].set_xlim(3,18)
		ax[1].set_xlim(3,18)
		ax[0].set_ylabel('Airmass')
		ax[1].set_ylabel(r'$v_{pl}$ [km/s]')
		ax[0].set_xlabel('UT time [hr]')
		ax[1].set_xlabel('UT time [hr]')
		ax[0].set_title('HST date: '+hstdates[j])
		plt.savefig('plots/'+hstdates[j]+'.png', bbox_inches='tight')
		plt.close()
		# plt.show()


