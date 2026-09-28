import numpy as np
import pytest
from astropy.table import Table

from hrccs_planner.cli import main
from hrccs_planner.contrast import dayside_contrast, dayside_temperature
from hrccs_planner.starlists import keck_lines, magellan_lines


def test_dayside_temperature():
    # no redistribution, zero albedo: T_day = T_eff sqrt(R*/a) (2/3)^1/4
    assert dayside_temperature(6000., 4., 0., 0.) == pytest.approx(3000. * (2 / 3.) ** 0.25)
    c_k, _ = dayside_contrast(6500., 1.5, 4., 1.5, 'K')
    c_h, _ = dayside_contrast(6500., 1.5, 4., 1.5, 'H')
    assert 0 < c_h < c_k < 0.01


def _simbad_row():
    return Table({'main_id': ['* tau Boo'], 'ra': [206.8156], 'dec': [17.4569], 'pmra': [-480.3],
                  'pmdec': [54.2], 'V': [4.5], 'R': [np.ma.masked]})


def test_keck_line_format():
    (line,) = keck_lines(['tau Boo'], _simbad_row())
    name, rest = line[:16], line[16:].split()
    assert name.strip() == 'tau Boo'
    assert rest[:7] == ['13', '47', '15.744', '+17', '27', '24.84', '2000']
    # pmra in seconds of time per year: mu_alpha* / (15000 cos dec)
    assert float(rest[7].split('=')[1]) == pytest.approx(-480.3 / 15000. / np.cos(np.radians(17.4569)), abs=1e-4)
    assert rest[-1] == 'Vmag=4.5'


def test_magellan_line_format():
    (line,) = magellan_lines(['tau Boo'], _simbad_row())
    f = line.split()
    assert f[:4] == ['001', 'tau_Boo', '13:47:15.74', '+17:27:24.84'] and line.endswith('# V=4.5')


def test_cli_init_and_windows(tmp_path, capsys):
    main(['init', str(tmp_path)])
    assert (tmp_path / 'targetlists' / 'examples_targetlist.csv').exists()
    main(['windows', 'examples', '--site', 'keck', '--start-date', '2027-07-10', '--end-date', '2027-07-11',
          '--workspace', str(tmp_path), '--csv', 'w.csv', '--split-targets'])
    assert (tmp_path / 'output' / 'best_windows' / 'w.csv').exists()
    assert (tmp_path / 'output' / 'best_windows' / 'w.txt').exists()
    assert list((tmp_path / 'output' / 'best_windows').glob('w_*.csv'))


def test_cli_missing_targetlist(tmp_path):
    with pytest.raises(SystemExit):
        main(['windows', 'nope', '--site', 'keck', '--start-date', '2027-07-10', '--end-date', '2027-07-10',
              '--workspace', str(tmp_path)])
