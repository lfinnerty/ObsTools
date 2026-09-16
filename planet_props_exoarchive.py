import argparse
import csv
import math
import re
import sys
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive


OUTPUT_COLUMNS = [
	'Name',
	'RA',
	'Dec',
	'Period (day)',
	'a (AU)',
	'T transit/inf conj (JD)',
	'Obs status',
	'Eccentric?',
	'Kp max (km/s)',
	'Transiting?',
	'Kmag (<9)',
	'e',
	'omega',
	'inc',
]

ARCHIVE_COLUMNS = (
	'hostname,pl_name,ra,dec,pl_orbper,pl_orbsmax,pl_tranmid,'
	'pl_orbeccen,pl_orblper,pl_orbincl,sy_kmag,tran_flag'
)

AU_PER_DAY_TO_KM_PER_SECOND = 1731.45683633


def _row_value(row, column):
	value = row[column]
	if value is None or getattr(value, 'mask', False):
		return None
	if hasattr(value, 'unit') and hasattr(value, 'value'):
		value = value.value
	try:
		value = value.item()
	except AttributeError:
		pass
	if isinstance(value, float) and not math.isfinite(value):
		return None
	return value


def _format_coordinate(ra, dec):
	coordinate = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg)
	formatted = coordinate.to_string('hmsdms', sep=':', precision=2, pad=True)
	return formatted.split()


def _format_number(value):
	if value is None:
		return ''
	return str(value)


def _make_output_row(row):
	ra, dec = _format_coordinate(_row_value(row, 'ra'), _row_value(row, 'dec'))
	period = _row_value(row, 'pl_orbper')
	semimajor_axis = _row_value(row, 'pl_orbsmax')
	eccentricity = _row_value(row, 'pl_orbeccen')

	kp_max = None
	if period is not None and semimajor_axis is not None:
		kp_max = (2 * math.pi * float(semimajor_axis) / float(period))
		### Note that the eccentricity correction is done later in the obs planning tool,
		### so it is commented out here
		# if eccentricity is not None and float(eccentricity) < 1:
		# 	kp_max /= math.sqrt(1 - float(eccentricity) ** 2)
		kp_max *= AU_PER_DAY_TO_KM_PER_SECOND

	transit_flag = _row_value(row, 'tran_flag')
	transiting = 'Y' if transit_flag == 1 else 'N'
	eccentric = 'Yes' if eccentricity is not None and float(eccentricity) > 0 else 'No'

	return {
		'Name': _row_value(row, 'pl_name'),
		'RA': ra,
		'Dec': dec,
		'Period (day)': _format_number(period),
		'a (AU)': _format_number(semimajor_axis),
		'T transit/inf conj (JD)': _format_number(_row_value(row, 'pl_tranmid')),
		'Obs status': 'None',
		'Eccentric?': eccentric,
		'Kp max (km/s)': _format_number(round(kp_max)) if kp_max is not None else '',
		'Transiting?': transiting,
		'Kmag (<9)': _format_number(_row_value(row, 'sy_kmag')),
		'e': _format_number(eccentricity),
		'omega': _format_number(_row_value(row, 'pl_orblper')),
		'inc': _format_number(_row_value(row, 'pl_orbincl')),
	}


def query_target_rows(system_names):
	conditions = []
	for name in system_names:
		name_variants = [name]
		spaced_name = re.sub(r'(?<=\d)([A-Z])(?= [a-z]$)', r' \1', name)
		if spaced_name != name:
			name_variants.append(spaced_name)
		for name_variant in name_variants:
			escaped_name = name_variant.replace("'", "''")
			conditions.append(
				"hostname = '{0}' OR pl_name = '{0}'".format(escaped_name)
			)

	rows = NasaExoplanetArchive.query_criteria(
		table='pscomppars',
		select=ARCHIVE_COLUMNS,
		where=' OR '.join(conditions),
	)
	def normalized_name(name):
		return re.sub(r'\s+', '', name).casefold()

	found_names = {
		normalized_name(str(value))
		for row in rows
		for value in (_row_value(row, 'hostname'), _row_value(row, 'pl_name'))
		if value is not None
	}
	missing_names = {
		name for name in system_names if normalized_name(name) not in found_names
	}
	if missing_names:
		raise ValueError(
			'No planets found for system(s): {}'.format(', '.join(sorted(missing_names)))
		)
	return rows


def output_path(output_name, system_names):
	if output_name is None:
		output_name = '_'.join(system_names)
	output_name = Path(output_name).name
	if output_name.lower().endswith('.csv'):
		output_name = output_name[:-4]
	if output_name.lower().endswith('_targetlist'):
		output_name = output_name[:-11]
	output_name = re.sub(r'[^A-Za-z0-9_.-]+', '_', output_name).strip('._')
	if not output_name:
		raise ValueError('output filename must contain at least one valid character')
	return Path('targetlists') / (output_name + '_targetlist.csv')


def main():
	parser = argparse.ArgumentParser(
		description='Generate a KPIC target list from NASA Exoplanet Archive systems.'
	)
	parser.add_argument(
		'system_names',
		nargs='+',
		help='one or more host system names, such as WASP-33 or HD 189733',
	)
	parser.add_argument(
		'-o', '--output',
		help='output filename; it is written under targetlists/ and gets a _targetlist.csv suffix',
	)
	args = parser.parse_args()

	try:
		rows = query_target_rows(args.system_names)
	except Exception as error:
		parser.error(str(error))

	try:
		output_file_path = output_path(args.output, args.system_names)
	except ValueError as error:
		parser.error(str(error))
	output_file_path.parent.mkdir(parents=True, exist_ok=True)
	with open(output_file_path, 'w', newline='') as output_file:
		writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
		writer.writeheader()
		for row in rows:
			writer.writerow(_make_output_row(row))

	print('Wrote {} planet(s) to {}'.format(len(rows), output_file_path))


if __name__ == '__main__':
	main()