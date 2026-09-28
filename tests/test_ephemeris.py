import numpy as np
import pytest

from hrccs_planner.ephemeris import (
    lambert_phase,
    orbital_phase,
    planet_rv_circular,
    planet_rv_eccentric,
    propagate_conjunction,
    solve_kepler,
)


def test_phase_and_circular_rv():
    assert orbital_phase(10.25, 10., 1.) == pytest.approx(0.25)
    assert orbital_phase(9.75, 10., 1.) == pytest.approx(0.75)
    # receding (positive) after transit, maximal at quadrature, zero at conjunctions
    assert planet_rv_circular(0.25, 150.) == pytest.approx(150.)
    assert planet_rv_circular([0., 0.5], 150.) == pytest.approx([0., 0.], abs=1e-9)


@pytest.mark.parametrize('e', [0., 0.1, 0.5])
def test_kepler(e):
    M = np.linspace(0, 2 * np.pi, 50)
    E = solve_kepler(M, e)
    assert np.allclose(E - e * np.sin(E), M, atol=1e-10)


def test_eccentric_reduces_to_circular():
    phase = np.linspace(0, 1, 30)
    # e = 0, omega = 90 deg: v = -Kp cos(2 pi phase + pi/2) = Kp sin(2 pi phase)
    assert np.allclose(planet_rv_eccentric(phase, 100., 0., 90.), planet_rv_circular(phase, 100.))


def test_lambert_phase():
    assert lambert_phase(0.5, 90.) == pytest.approx(1.)       # full phase at secondary eclipse
    assert lambert_phase(0., 90.) == pytest.approx(0.)        # new phase at transit
    assert lambert_phase(0.25, 90.) == pytest.approx(1 / np.pi)
    # at i = 30 deg the phase angle only spans 60-120 deg
    a_min, a_max = np.radians(60.), np.radians(120.)
    f = lambda a: (np.sin(a) + (np.pi - a) * np.cos(a)) / np.pi
    assert lambert_phase(0.5, 30.) == pytest.approx(f(a_min))
    assert lambert_phase(0., 30.) == pytest.approx(f(a_max))


def test_propagate_conjunction():
    n, tn, st, sp = propagate_conjunction(1000., 0.01, 2., 0.001, 1000. + 2. * 300 + 0.3)
    assert n == 300 and tn == pytest.approx(1600.)
    assert st == pytest.approx(np.hypot(0.01, 0.3)) and sp == pytest.approx(st / 2.)
