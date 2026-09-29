"""hrccs_planner: observation planning for high-resolution cross-correlation
spectroscopy (HRCCS) of exoplanet atmospheres."""
__version__ = '0.1.0'

from .ephemeris import (  # noqa: F401
    eccentric_orbit,
    eclipse_phase,
    lambert_phase,
    orbital_phase,
    planet_rv_circular,
    planet_rv_eccentric,
    propagate_conjunction,
)
from .sites import SITES, Site, get_site  # noqa: F401
from .targets import Target, read_targets, write_targets  # noqa: F401
