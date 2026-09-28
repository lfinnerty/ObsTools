import os

import pytest

from hrccs_planner.targets import read_targets
from hrccs_planner.workspace import EXAMPLES

DATA = os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture(scope='session')
def examples():
    return read_targets(os.path.join(EXAMPLES, 'examples_targetlist.csv'))
