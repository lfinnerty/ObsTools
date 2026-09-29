import dataclasses

import numpy as np
import pytest

from hrccs_planner.cli import main
from hrccs_planner.ephemeris import propagate_conjunction
from hrccs_planner.phase_sensitivity import (
    DELTA_GRID,
    Scenario,
    expected_ratio,
    phase_sensitivity,
    sigma_phase,
    snr_versus_shift,
)
from hrccs_planner.sites import get_site
from hrccs_planner.visibility import dates_between
from hrccs_planner.windows import WindowOptions, rank_windows

DATES = dates_between('2027-07-10', '2027-07-19')
OPTS = WindowOptions(duration=2., normalize_per_target=True, snr_min=0., sort='snr')


@pytest.fixture(scope='module')
def kelt9(examples):
    target = next(t for t in examples if t.name == 'KELT-9 b')
    return dataclasses.replace(target, extra={**target.extra, 'sig 0 h': '0', 'sig 1 h': '1', 'sig 3 h': '3'})


def test_expected_ratio_limits():
    ratio = 1 - 5 * DELTA_GRID ** 2
    assert expected_ratio(ratio, 0.) == (1., 1.)
    assert all(np.isnan(expected_ratio(ratio, np.nan)))
    means = [expected_ratio(ratio, s)[0] for s in (0.005, 0.01, 0.02, 0.04)]
    assert np.all(np.diff(means) < 0)
    # Second-order loss: 1 - 5 sigma^2 for a quadratic ratio.
    assert means[2] == pytest.approx(1 - 5 * 0.02 ** 2, rel=1e-4)


def test_snr_is_unchanged_without_a_shift_and_falls_off_eclipse(kelt9):
    site = get_site('keck')
    rows, _ = rank_windows([kelt9], site, DATES, OPTS, log=lambda s: None)
    best = max(rows, key=lambda r: r['snr'])
    ratio = snr_versus_shift(kelt9, best, site, OPTS)
    assert ratio[np.argmin(np.abs(DELTA_GRID))] == pytest.approx(1.)
    assert ratio.max() < 1.1
    assert ratio[np.abs(DELTA_GRID) > 0.15].max() < 0.9


def test_sigma_phase_sources(kelt9):
    jd = kelt9.t0 + 500.3
    assert sigma_phase(kelt9, Scenario('three', 'sig 3 h'), jd) == pytest.approx(3 / 24 / kelt9.period)
    assert np.isnan(sigma_phase(kelt9, Scenario('missing', 'no such column'), jd))
    with_errors = dataclasses.replace(kelt9, extra={'T err (d)': '0.001', 'P err (d)': '1e-5'})
    expected = propagate_conjunction(kelt9.t0, 0.001, kelt9.period, 1e-5, jd)[3]
    assert sigma_phase(with_errors, Scenario('ephemeris'), jd) == pytest.approx(expected)
    assert sigma_phase(kelt9, Scenario('ephemeris'), jd) == 0.   # no error columns


def test_better_ephemeris_recovers_snr(kelt9):
    scenarios = [Scenario(c, c) for c in ('sig 3 h', 'sig 1 h', 'sig 0 h')]
    summary, windows = phase_sensitivity([kelt9], get_site('keck'), DATES, OPTS, scenarios, top=3, log=lambda s: None)
    means = {s['scenario']: s['mean_ratio'] for s in summary}
    assert means['sig 0 h'] == 1.
    assert means['sig 3 h'] < means['sig 1 h'] < 1.
    assert len(windows) == 3 * len(scenarios)
    assert {round(s['sigma_t_h'], 6) for s in summary} == {0., 1., 3.}


def test_cli_phase_sensitivity(tmp_path, capsys):
    targets = tmp_path / 'targets.csv'
    targets.write_text(
        'Name,RA,Dec,Period (day),T transit/inf conj (JD),Kp max (km/s),Eccentric?,sig now (h),sig new (h)\n'
        'KELT-9 b,20:31:26.38,+39:56:20.10,1.4811235,2457095.68572,247,No,2,0.5\n')
    main(['phase-sensitivity', str(targets), '--site', 'keck', '--start', '2027-07-10', '--end', '2027-07-19',
          '--duration', '2', '--sigma-column', 'sig now (h)', '--sigma-column', 'sig new (h)',
          '--top-windows', '2', '--csv', str(tmp_path / 'ps.csv')])
    out = capsys.readouterr().out
    assert 'KELT-9 b' in out and "relative to 'sig now (h)'" in out
    assert (tmp_path / 'ps.csv').exists() and (tmp_path / 'ps_windows.csv').exists()
    lines = (tmp_path / 'ps.csv').read_text().splitlines()
    assert lines[0].startswith('target,scenario') and len(lines) == 3
