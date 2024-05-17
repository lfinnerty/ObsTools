from astroquery.simbad import Simbad
import numpy as np
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
	datestr = '2024-05-19'
	namestr = 'KPIC'
	tgtnames = ['HIP 62944',
				'HIP 81497',
				'HIP 61960',
				'ups Her',
				'WASP-14',
				'HIP 56147',
				'HIP 55507'
			   ]
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
