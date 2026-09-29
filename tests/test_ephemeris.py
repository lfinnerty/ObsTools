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


@pytest.mark.parametrize('omega', [0., 90., 200., -85.])
def test_eccentric_reduces_to_circular(omega):
    from hrccs_planner.ephemeris import eccentric_orbit, eclipse_phase
    phase = np.linspace(0, 1, 30)
    assert np.allclose(planet_rv_eccentric(phase, 100., 0., omega), planet_rv_circular(phase, 100.))
    uphase, rscale, _ = eccentric_orbit(phase, 0., omega)
    assert np.allclose(uphase, phase) and np.allclose(rscale, 1.)
    assert eclipse_phase(0., omega) == pytest.approx(0.5)


@pytest.mark.parametrize('e,omega', [(0.3, 60.), (0.6, 200.), (0.93, -58.9)])
def test_eccentric_orbit_against_kepler(e, omega):
    from hrccs_planner.ephemeris import eccentric_orbit, eclipse_phase
    # independent orbit: planet position from Kepler's equation (Newton), T0 = transit
    w = np.radians(omega)
    f0 = np.pi / 2 - w
    E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f0 / 2))
    tp = -(E0 - e * np.sin(E0)) / (2 * np.pi)          # periastron, in periods from T0
    phase = np.linspace(-0.2, 1.2, 14001)
    M = 2 * np.pi * (phase - tp)
    E = M + e * np.sin(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    f = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = (1 - e ** 2) / (1 + e * np.cos(f))
    z = r * np.sin(w + np.pi + f)                        # line of sight, away from us
    x = r * np.cos(w + np.pi + f)                        # along the motion at transit
    uphase, rscale, rv = eccentric_orbit(phase, e, omega)
    # velocity in units of 2 pi a / P: the semi-amplitude is 1/sqrt(1-e^2)
    assert np.allclose(rv, np.gradient(z, phase) / (2 * np.pi), atol=2e-3 * np.max(np.abs(rv)))
    # (solve_kepler's fixed-point iteration is good to ~1e-4 at e = 0.93)
    assert np.allclose(rscale * np.sin(2 * np.pi * uphase), x, atol=5e-4)
    assert np.allclose(-rscale * np.cos(2 * np.pi * uphase), z, atol=5e-4)
    # at the eclipse time the planet is behind the star (x = 0, z > 0)
    uph_e, rs_e, _ = eccentric_orbit(eclipse_phase(e, omega), e, omega)
    assert (uph_e % 1) == pytest.approx(0.5, abs=1e-4)


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
