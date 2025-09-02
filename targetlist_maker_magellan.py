### Turning Lennart's juptyer notebook into a script

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def make_line(refnumber, targetname, ra, dec, eq, RApm, decPM, rot_angle='0.0', 
	    rot_mode='GRV', gs1_ra='0', gs1_dec='0', gs1_eq='2000.0', gs2_ra='0', gs2_dec='0', gs2_eq='2000.0', epoch='2017.5'):
	line = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
	    refnumber, targetname, ra, dec, eq, RApm, decPM, rot_angle, 
	    rot_mode, gs1_ra, gs1_dec, gs1_eq, gs2_ra, gs2_dec, gs2_eq, epoch)
	return line

if __name__ == '__main__':
	datestr = '2025-09-29'
	namestr = 'winered'
	tgtnames = ['WASP-74',
				'WASP-94B',
				'WASP-18',
				'WASP-77A',
				'WASP-121']
	

	hdr = "\n###".join(("""# RA Dec equinox RApm Decpm offset rot RA_probe1 Dec_probe1 equinox RA_probe2 Dec_probe2 equinox pm_epoch""",
              """name hh:mm:ss.s sdd:mm:ss yyyy.0 s.ss s.ss angle mode hh:mm:ss.s sdd:mm:ss yyyy.0 hh:mm:ss.s sdd:mm:ss yyyy.0yyyy.0\n"""))

	# colwidth = 16

	sbad = Simbad()
	sbad.add_votable_fields('pmra')
	sbad.add_votable_fields('pmdec')
	sbad.add_votable_fields('fluxdata(V)')
	result_table = sbad.query_objects(tgtnames)

	outfile = open('starlists/targetlist_magellan_'+namestr+'_'+datestr+'.txt', 'w+')
	for i, entry in enumerate(result_table):
		if entry['main_id'] == '':
			print('Error, invalid query on', tgtnames[i])
			raise RuntimeError
		
		if len(tgtnames[i]) > 16:
			tgtnames[i] = tgtnames[i][:15]
		
		coord = SkyCoord(ra=entry['ra']*u.deg,dec=entry['dec']*u.deg)
		cstring = coord.to_string('hmsdms')
		
		ra, dec = cstring.split(' ')
		ra = ra.replace('h',':')
		ra = ra.replace('m',':')
		ra = ra[:11]

		dec = dec.replace('d',':')
		dec = dec.replace('m',':')
		dec = dec[:12]



		pmra = str(np.round(entry['pmra']/1e3 / (15.*np.cos(np.pi/180. * float(dec.split(':')[0]))),decimals=3))
		pmdec =str(np.round(entry['pmdec']/1e3,decimals=3))

		result_string = make_line(str(i+1).zfill(3),tgtnames[i],ra,dec,'2000.0',pmra,pmdec)

		if entry['V'] != '--':
			vmag = ' # V='+str(np.round(entry['V'],1))
			result_string+=vmag
		print(result_string)
		outfile.write(result_string+'\n')
	outfile.close()
