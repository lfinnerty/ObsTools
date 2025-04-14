from astroquery.simbad import Simbad
import numpy as np
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
	datestr = '2025-01-10'
	namestr = 'KPIC'
	tgtnames = ['HD 35067',
				'HIP 19388',
				'HIP 81497',
				'HIP 95771',
				'HIP 62944',
				'HIP 110882',
				'79 Cyg',
				'HIP 61960',
				'HD 130163',
				'* z cen',
				'ups Her',
				'HIP 82714',
				'HAT-P-22',
				'HD 88133',
				'GJ 436',
				'HD 49674',
				'55 Cnc',
				'WASP-33',
				'KELT-2A',
				'KELT-7',
				'WASP-127',
				'HD 102195',
				'BD-10 3166',
				'RX J0342.5+1216',
				'HIP16322',
				'DH Tau',
				'HIP18717']
	colwidth = 16

	sbad = Simbad()
	sbad.add_votable_fields('pm')
	# sbad.add_votable_fields('fluxdata(K)')
	sbad.add_votable_fields('fluxdata(R)')
	sbad.add_votable_fields('fluxdata(V)')
	result_table = sbad.query_objects(tgtnames)
	# print(result_table)
	if len(result_table.errors) != 0:
		print('Errors encountered!')
		for error in result_table.errors:
			print(error)
		raise RuntimeError('Invalid query!')

	outfile = open('starlists/starlist_'+namestr+'_'+datestr+'.txt', 'w+')
	for i, entry in enumerate(result_table):
		# print(entry.keys())
		# print(entry['FLUX_K'], entry['FLUX_R'], entry['FLUX_V'])
		# result_string = entry['MAIN_ID'].ljust(colwidth, ' ')
		if len(tgtnames[i]) > 16:
			tgtnames[i] = tgtnames[i][:15]
		result_string = tgtnames[i].ljust(colwidth, ' ')
		result_string+=entry['RA'].ljust(colwidth, ' ')
		result_string+=entry['DEC'].ljust(colwidth, ' ')
		result_string+= '2000 '

		pmra = 'pmra='+str(np.round(entry['PMRA']/1e3 / 15.,decimals=4))
		pmdec = 'pmdec='+str(np.round(entry['PMDEC']/1e3,decimals=4))

		result_string+= pmra.ljust(colwidth, ' ') + pmdec.ljust(colwidth, ' ' )
		# if entry['FLUX_K'] != '--':
		# 	kmag = 'Kmag='+str(entry['FLUX_K'])
		# 	result_string+=kmag.ljust(colwidth, ' ')
		if entry['FLUX_R'] != '--':
			rmag = 'Rmag='+str(entry['FLUX_R'])
			result_string+=rmag.ljust(colwidth, ' ')
		elif entry['FLUX_V'] != '--':
			vmag = 'Vmag='+str(entry['FLUX_V'])
			result_string+=vmag.ljust(colwidth, ' ')
		print(result_string)
		outfile.write(result_string+'\n')
	outfile.close()
