from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
	datestr = '2025-06-26'
	namestr = 'KPIC'
	# tgtnames = ['HD 35067',
	# 			'HIP 19388',
	# 			'HIP 81497',
	# 			'HIP 95771',
	# 			'HIP 62944',
	# 			'HIP 110882',
	# 			'79 Cyg',
	# 			'HIP 61960',
	# 			'HD 130163',
	# 			'* z cen',
	# 			'ups Her',
	# 			'HIP 82714',
	# 			'HAT-P-22',
	# 			'HD 88133',
	# 			'GJ 436',
	# 			'HD 49674',
	# 			'55 Cnc',
	# 			'WASP-33',
	# 			'KELT-2A',
	# 			'KELT-7',
	# 			'WASP-127',
	# 			'HD 102195',
	# 			'BD-10 3166',
	# 			'RX J0342.5+1216',
	# 			'HIP16322',
	# 			'DH Tau',
	# 			'HIP18717']
	tgtnames = ['HIP 81497',
				'HIP 95771',
				'TOI 1518',
				'MASCARA-1']
	colwidth = 16

	sbad = Simbad()
	sbad.add_votable_fields('pmra')
	sbad.add_votable_fields('pmdec')
	sbad.add_votable_fields('fluxdata(R)')
	sbad.add_votable_fields('fluxdata(V)')
	result_table = sbad.query_objects(tgtnames)

	outfile = open('starlists/starlist_'+namestr+'_'+datestr+'.txt', 'w+')
	for i, entry in enumerate(result_table):
		if entry['main_id'] == '':
			print('Error, invalid query on', tgtnames[i])
			raise RuntimeError
		
		if len(tgtnames[i]) > 16:
			tgtnames[i] = tgtnames[i][:15]
		
		coord = SkyCoord(ra=entry['ra']*u.deg,dec=entry['dec']*u.deg)
		cstring = coord.to_string('hmsdms')
		
		ra, dec = cstring.split(' ')
		ra = ra.replace('h',' ')
		ra = ra.replace('m',' ')
		ra = ra[:12]

		dec = dec.replace('d',' ')
		dec = dec.replace('m',' ')
		dec = dec[:12]

		result_string = tgtnames[i].ljust(colwidth, ' ')
		result_string+=ra.ljust(colwidth, ' ')
		result_string+=dec.ljust(colwidth,' ')
		result_string+= '2000 '

		pmra = 'pmra='+str(np.round(entry['pmra']/1e3 / 15.,decimals=4))
		pmdec = 'pmdec='+str(np.round(entry['pmdec']/1e3,decimals=4))

		result_string+= pmra.ljust(colwidth, ' ') + pmdec.ljust(colwidth, ' ' )

		if entry['R'] != '--':
			rmag = 'Rmag='+str(np.round(entry['R'],2))
			result_string+=rmag.ljust(colwidth, ' ')
		elif entry['V'] != '--':
			vmag = 'Vmag='+str(np.round(entry['V'],2))
			result_string+=vmag.ljust(colwidth, ' ')
		print(result_string)
		outfile.write(result_string+'\n')
	outfile.close()
