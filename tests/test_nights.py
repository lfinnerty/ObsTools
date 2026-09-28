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
