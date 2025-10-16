import numpy as np 
import matplotlib.pyplot as plt 
import sys
import configparser
import ast
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

h = 6.626e-27
c = 3e10
kb = 1.38e-16

def flux_ratio(waves, t1,t2,r1,r2):
	waves = waves/1e4
	b1 = 1/(np.exp(h*c/(waves*kb*t1)) - 1)
	b2 = 1/(np.exp(h*c/(waves*kb*t2)) - 1)

	return b1/b2 * r1**2/r2**2

if __name__ == '__main__':
	inst = sys.argv[1]
	objects = np.genfromtxt(inst+'_targetlist.csv', delimiter=',', skip_header=1, dtype=None, encoding=None)
	names = ['WASP-76', 'WASP-77', 'WASP-121', 'WASP-122', 'WASP-87', 'WASP-103', 'MASCARA-1', 'WASP-19', 'WASP-82', 'MASCARA-4', 'WASP-189', 'WASP-74']


	# snr_st_5 = np.array([241, 208,138,124,124,57,338,65,188,338,750,241])
	snr_st_5 = np.array([447,388,265,241,241,124,620,138,353,620,1300,447])
	tint_min = 0.9*np.array([300,300,300,300,300,300,300,300,300,300,300,300])
	snrs = snr_st_5*np.sqrt(tint_min/5)
	# print(snrs)

	nlines = np.array([2000,2000,400,2000,400,400,400,2000,400,400,70,2000])	
	eta = 0.55
	
	for i, name in enumerate(names):
		for j in range(len(objects)):
			if objects[j][0] == name:
				obj = objects[j]

		tst = float(obj[17])
		rst = float(obj[18])
		tpl = float(obj[16])
		rpl = float(obj[15]*0.10045)
		fpfs = flux_ratio(1.2,tpl,tst,rpl,rst)

		print(name, tpl, nlines[i], eta*fpfs*snrs[i]*np.sqrt(nlines[i]))


