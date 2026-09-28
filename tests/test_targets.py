import math

import pytest

from hrccs_planner.targets import STANDARD_HEADER, Target, TargetListError, read_targets, write_targets

LEGACY = '''Name,RA,Dec,Period (day),a (AU),T transit/inf conj (JD),KPIC obs status,Eccentric?,Kp max (km/s),Transiting?,Kmag (<9),e,omega,inc,Mpl,People interested
"HD 1, the star",09:22:37.67,+50:36:13.6,111.436765,0.46,2458888.07466,None,Yes,123.8,Yes,7,0.93183,-58.887,90,4.16,Someone
WASP-33 b,02:26:51.06,+37:33:01.7,1.2198675,0.0259,2458804.832638889,None,No,230.9,Y,,0,0,
'''


def test_read_legacy_header(tmp_path):
    p = tmp_path / 'x_targetlist.csv'
    p.write_text(LEGACY)
    a, b = read_targets(p)
    assert a.name == 'HD 1, the star' and a.eccentric and a.transiting and a.e == pytest.approx(0.93183)
    assert a.status == 'None' and a.extra == {'Mpl': '4.16', 'People interested': 'Someone'}
    assert not b.eccentric and math.isnan(b.kmag) and b.inc is None
    assert b.coord.ra.deg == pytest.approx(36.71275)
    assert b.host == 'WASP-33'


def test_missing_column(tmp_path):
    p = tmp_path / 'bad.csv'
    p.write_text('Name,RA,Dec,Period (day)\nX,1,2,3\n')
    with pytest.raises(TargetListError, match='Kp max'):
        read_targets(p)


def test_bad_value_reports_line(tmp_path):
    p = tmp_path / 'bad.csv'
    p.write_text(LEGACY.replace('1.2198675', 'one'))
    with pytest.raises(TargetListError, match='line 3'):
        read_targets(p)


def test_roundtrip(tmp_path):
    t = Target('P b', '10:00:00', '+10:00:00', 2.5, 2460000.5, 150., kmag=8.1, extra={'Notes': 'n'})
    p = tmp_path / 'rt.csv'
    write_targets(p, [t])
    assert p.read_text().splitlines()[0].split(',') == STANDARD_HEADER + ['Notes']
    (u,) = read_targets(p)
    assert (u.name, u.period, u.kp_max, u.kmag, u.eccentric, u.extra) == ('P b', 2.5, 150., 8.1, False, {'Notes': 'n'})


def test_examples_load(examples):
    assert len(examples) == 9
    assert {t.name for t in examples} >= {'WASP-121 b', 'WASP-33 b', 'KELT-9 b', 'HD 80606 b'}
