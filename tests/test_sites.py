import numpy as np
import pytest

from hrccs_planner.sites import get_site, keck2_min_altitude


@pytest.mark.parametrize('name,key', [('KECK', 'maunakea'), ('Gemini-N', 'maunakea'), ('magellan', 'lco'),
                                      ('VLT', 'paranal'), ('gemini_south', 'gemini-s'), ('keck2', 'keck2')])
def test_aliases(name, key):
    assert get_site(name).key == key


def test_night_grid_covers_night():
    # the grid for local date D at Maunakea runs from D 21:00 UT to D+1 21:00 UT (11:00 HST to 11:00 HST)
    times, hours = get_site('keck').night_times('2027-07-10', 1441)
    assert hours[0] == -3 and hours[-1] == 21 and len(times) == 1441
    assert times[0].isot.startswith('2027-07-10T20:59:59')
    _, hours = get_site('lco').night_times('2027-01-05', 100)
    assert hours[0] == -6


def test_keck2_limit():
    assert np.all(keck2_min_altitude([0., 200., 300., 340.]) == [18., 36.8, 36.8, 18.])


def test_unknown_site():
    with pytest.raises(ValueError, match='unknown site'):
        get_site('not an observatory')
