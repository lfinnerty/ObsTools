"""Orbital phase, planet radial velocity, phase curve and ephemeris propagation.

Conventions:
    - phase 0 is the transit / inferior conjunction T0, phase 0.5 the secondary
      eclipse / superior conjunction (circular orbits);
    - the planet's line-of-sight velocity (relative to the star) is
      v = Kp sin(2 pi phase) for a circular orbit, positive (receding) just after
      transit, where Kp is the observed semi-amplitude Kp_max sin(i);
    - Kp_max = 2 pi a / P is the total orbital velocity.
"""
import numpy as np


def orbital_phase(jd, t0, period):
    """Orbital phase in [0, 1) at times `jd` (days, same time system as t0)."""
    return ((np.asarray(jd) - t0) / period) % 1


def planet_rv_circular(phase, kp):
    """Planet line-of-sight velocity [km/s] on a circular orbit."""
    return kp * np.sin(2 * np.pi * np.asarray(phase))


def solve_kepler(mean_anomaly, e, order=100):
    """Eccentric anomaly E from the mean anomaly M (fixed-point iteration of
    E = M + e sin E, `order` times)."""
    E = np.asarray(mean_anomaly, dtype=float)
    for _ in range(order):
        E = mean_anomaly + e * np.sin(E)
    return E


def planet_rv_eccentric(phase, kp, e, omega_deg):
    """Planet line-of-sight velocity [km/s] on an eccentric orbit.

    `phase` is (t - T0)/P, i.e. the mean anomaly over 2 pi measured from T0, and
    omega is the argument of periastron [deg]; v = -Kp (cos(f + omega) + e cos(omega)),
    with the true anomaly f from Kepler's equation.
    """
    omega = np.radians(omega_deg)
    E = solve_kepler(2 * np.pi * np.asarray(phase), e)
    beta = e / (1 + np.sqrt(1 - e ** 2))
    f = E + 2 * np.arctan(beta * np.sin(E) / (1 - beta * np.cos(E)))
    return -kp * (np.cos(f + omega) + e * np.cos(omega))


def lambert_phase(phase, inclination_deg=90.):
    """Lambertian phase curve, normalized to 1 at full phase.

    f = (sin a + (pi - a) cos a) / pi with the phase angle a given by
    cos a = -sin(i) cos(2 pi phase); for an edge-on orbit a = pi |1 - 2 phase|
    (a = 0 at secondary eclipse).
    """
    sini = np.sin(np.radians(inclination_deg))
    alpha = np.arccos(np.clip(-sini * np.cos(2 * np.pi * np.asarray(phase)), -1, 1))
    return (np.sin(alpha) + (np.pi - alpha) * np.cos(alpha)) / np.pi


def propagate_conjunction(t0, t0_err, period, period_err, t):
    """Nearest conjunction to time `t` and its uncertainty.

    Returns (n, t_n, sigma_t, sigma_phase) with t_n = t0 + n P and
    sigma_t = sqrt(t0_err^2 + (n period_err)^2) (uncorrelated T0 and P errors).
    """
    n = np.round((np.asarray(t) - t0) / period)
    sigma_t = np.sqrt(t0_err ** 2 + (n * period_err) ** 2)
    return n, t0 + n * period, sigma_t, sigma_t / period
