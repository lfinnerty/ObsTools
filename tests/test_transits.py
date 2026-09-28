import numpy as np
import pytest

from hrccs_planner.sites import get_site
from hrccs_planner.targets import Target
from hrccs_planner.transits import TransitOptions, choose_baseline, rank_transits
from hrccs_planner.visibility import dates_between, night_sky

HOURS = np.linspace(0., 10., 601)  # 1-minute grid


def _ok(a, b):
    return (HOURS >= a) & (HOURS <= b)


def test_split_baseline_needs_both_sides():
    ok = _ok(2., 8.)
    assert choose_baseline(ok, HOURS, 4., 5.5, 2., 'split') == pytest.approx((1., 1.))
    # usable time ends 0.5 h after egress: a split 1 h + 1 h baseline no longer fits
    ok = _ok(2., 6.)
    assert choose_baseline(ok, HOURS, 4., 5.5, 2., 'split') is None


def test_any_baseline_most_even_that_fits():
    ok = _ok(2., 6.)
    pre, post = choose_baseline(ok, HOURS, 4., 5.5, 2., 'any')
    assert pre + post == pytest.approx(2.)
    assert post == pytest.approx(0.5, abs=0.02) and pre == pytest.approx(1.5, abs=0.02)
    # all on one side when the transit sits at the end of the usable time
    pre, post = choose_baseline(_ok(1., 5.5), HOURS, 4., 5.5, 2., 'any')
    assert (pre, post) == pytest.approx((2., 0.))
    assert choose_baseline(_ok(3., 5.5), HOURS, 4., 5.5, 2., 'any') is None


def test_windows_meet_the_constraints(examples):
    site = get_site('keck')
    opts = TransitOptions(duration=3., baseline=2., baseline_mode='any', snr_min=0., airmass_k=0.05)
    rows, _ = rank_transits(examples, site, dates_between('2027-07-10', '2027-07-19'), opts, log=lambda s: None)
    assert rows
    by_name = {t.name: t for t in examples}
    for r in rows:
        assert r['pre'] + r['post'] == pytest.approx(2.)
        assert r['ut1'] - r['ut0'] == pytest.approx(5.)
        night = night_sky(site, r['date'])
        alt = by_name[r['name']].coord.transform_to(night.frame).alt.deg
        win = (night.hours >= r['ut0']) & (night.hours <= r['ut1'])
        assert np.all(alt[win] > 30.) and np.all(night.sun_alt[win] < -18.)
        assert r['X_max'] < 2.
        # mid-transit agrees with T0 + n P
        t = by_name[r['name']]
        jd_mid = night.times.jd[0] + (r['tc'] - night.hours[0]) / 24.
        n = (jd_mid - t.t0) / t.period
        assert abs(n - round(n)) * t.period * 24 < 1e-6


def test_non_transiting_and_missing_duration_are_skipped():
    t = Target('X b', '10:00:00', '+20:00:00', 2., 2460000., 150., transiting=False)
    msgs = []
    rows, _ = rank_transits([t], get_site('keck'), ['2027-03-01'], TransitOptions(duration=2.), log=msgs.append)
    assert not rows and 'not transiting' in msgs[0]
    t.transiting = True
    rows, _ = rank_transits([t], get_site('keck'), ['2027-03-01'], TransitOptions(), log=msgs.append)
    assert 'no transit duration' in msgs[-1]
