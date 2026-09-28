import json
import os

import pytest
from conftest import DATA

from hrccs_planner.nights import NightOptions, observable_targets
from hrccs_planner.sites import get_site
from hrccs_planner.visibility import dates_between

REF = json.load(open(os.path.join(DATA, 'nights_reference.json')))


@pytest.mark.parametrize('site,dates,kw', [
    ('keck', ('2027-07-10', '2027-07-16'), {}),
    ('lco', ('2027-01-05', '2027-01-09'), dict(require_conjunction=True, inclination=60.)),
])
def test_regression_against_legacy(examples, site, dates, kw):
    """Same targets on the same nights as the original observing_planner.py."""
    key = 'keck' if site == 'keck' else 'lco'
    got = {}
    for date in dates_between(*dates):
        _, shown = observable_targets(examples, get_site(site), date, NightOptions(**kw), log=lambda s: None)
        if shown:
            got[date] = sorted(t.name for t, *_ in shown)
    assert got == REF[key]


def test_parse_window_and_shading(tmp_path):
    from hrccs_planner.nights import parse_window, read_shading
    assert parse_window('07:49-14:33') == pytest.approx((7 + 49 / 60., 14.55))
    assert parse_window('22:30-02:00') == pytest.approx((-1.5, 2.))      # across 00 UT
    p = tmp_path / 't.csv'
    p.write_text('rank,hst_date,target,ut_mid,t14_hrs,ut_start,ut_end\n1,2027-02-24,X b,8.9,1.5,7.36,10.45\n')
    assert read_shading([p]) == {('2027-02-24', 'X b'): [(7.36, 10.45, pytest.approx(8.15), pytest.approx(9.65))]}


def test_compact_plot_with_shading(tmp_path, examples):
    import matplotlib
    matplotlib.use('Agg')
    from hrccs_planner.nights import plan_nights
    wasp33 = [t for t in examples if t.name == 'WASP-33 b']
    shading = {('2027-07-10', 'WASP-33 b'): [(12., 15., None, None)]}
    saved = plan_nights(wasp33, get_site('keck'), ['2027-07-10'], str(tmp_path), log=lambda s: None,
                        shading=shading, compact=True, fmt='pdf')
    assert len(saved) == 1 and saved[0].endswith('_compact.pdf') and os.path.getsize(saved[0]) > 0


def test_shaded_targets_are_always_shown(examples):
    """A transit window is on the night side of the orbit, which the dayside criterion rejects."""
    t = [x for x in examples if x.name == 'KELT-9 b']
    date = '2027-07-10'
    _, shown = observable_targets(t, get_site('keck'), date, NightOptions(delta_v_min=1e4), log=lambda s: None)
    assert not shown
    _, shown = observable_targets(t, get_site('keck'), date, NightOptions(delta_v_min=1e4), log=lambda s: None,
                                  always={'KELT-9 b'})
    assert [x.name for x, *_ in shown] == ['KELT-9 b']


def test_star_label():
    from hrccs_planner.nights import star_label
    assert star_label('TOI-1408') == 'TOI-1408A' and star_label('55 Cnc') == '55 Cnc A'
