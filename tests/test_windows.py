import os

import numpy as np
import pytest
from conftest import DATA

from hrccs_planner.sites import get_site
from hrccs_planner.visibility import dates_between
from hrccs_planner.windows import WindowOptions, best_subwindow, rank_windows, write_windows


def test_best_subwindow_picks_the_best_stretch():
    usable = np.zeros(100, bool)
    usable[10:90] = True
    w2 = np.ones(100)
    w2[60:70] = 10.
    vpl = np.linspace(100, -100, 100)
    sub = best_subwindow(usable, w2, vpl, 20)
    assert sub.sum() == 20 and sub[60:70].all()
    assert np.array_equal(best_subwindow(usable, w2, vpl, 200), usable)


# (reference CSV from the original best_window.py, site, dates, options)
CASES = [
    ('windows_keck_default.csv', 'KECK', ('2027-07-10', '2027-07-19'), dict(duration=2., snr_min=0.)),
    ('windows_gemn_cap_airmass.csv', 'GEMINI-N', ('2027-07-10', '2027-07-19'),
     dict(inclination=30., weight_kmag=True, max_hours=4., airmass_k=0.05, seeing_exponent=0.6, snr_min=0.)),
    ('windows_lco_pertarget.csv', 'LCO', ('2027-01-05', '2027-01-14'),
     dict(duration=2.5, alt_min=35., normalize_per_target=True, snr_min=0.5, top=40)),
]


@pytest.mark.parametrize('ref,site,dates,kw', CASES)
def test_regression_against_legacy(tmp_path, examples, ref, site, dates, kw):
    rows, _ = rank_windows(examples, get_site(site), dates_between(*dates), WindowOptions(**kw), log=lambda s: None)
    out = tmp_path / 'w.csv'
    write_windows(str(out), rows, '', command='')
    assert out.read_text() == open(os.path.join(DATA, ref)).read()


def test_eccentric_targets_are_skipped(examples):
    msgs = []
    rank_windows(examples, get_site('keck'), ['2027-07-10'], log=msgs.append)
    assert any('HD 80606 b' in m and 'eccentric' in m for m in msgs)
