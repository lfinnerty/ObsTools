"""Orbital phase, planet radial velocity, phase curve and ephemeris propagation.

Conventions:
    - phase 0 is the transit / inferior conjunction T0, phase 0.5 the secondary
      eclipse / superior conjunction (circular orbits). For an eccentric orbit T0 is
      still the transit time; the orbital phase u / 2 pi, with u = f + omega - pi/2
      the orbital angle from transit, is 0.5 at secondary eclipse and replaces the
      time phase for the dayside, eclipse and phase-curve geometry;
    - the planet's line-of-sight velocity (relative to the star) is
      v = Kp sin(2 pi phase) for a circular orbit, positive (receding) just after
      transit, where Kp is the observed semi-amplitude Kp_max sin(i);
    - Kp_max = 2 pi a / P is the total orbital velocity; an eccentric orbit's
      semi-amplitude is Kp / sqrt(1 - e^2);
    - omega is the star's argument of periastron (the RV convention), in degrees.

The eccentric-orbit functions match KPIC-Synthgen and KPIC-XCORR-Rebase.
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


def transit_mean_anomaly(e, omega_deg):
    """Mean anomaly [rad, in [0, 2 pi)] at transit, where the true anomaly is
    f = pi/2 - omega."""
    f = np.pi / 2 - np.radians(omega_deg)
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f / 2), np.sqrt(1 + e) * np.cos(f / 2))
    return (E - e * np.sin(E)) % (2 * np.pi)


def eclipse_phase(e, omega_deg):
    """Time of secondary eclipse (f = 3 pi/2 - omega) after transit, as a fraction of
    the period in [0, 1); 0.5 for a circular orbit."""
    return ((transit_mean_anomaly(e, omega_deg + 180.) - transit_mean_anomaly(e, omega_deg))
            % (2 * np.pi)) / (2 * np.pi)


def eccentric_orbit(phase, e, omega_deg):
    """Orbital phase, separation and RV factor on an eccentric orbit.

    `phase` is the time since transit T0 in periods, (t - T0)/P.

    Returns (uphase, rscale, rv): the orbital phase u / 2 pi (u = f + omega - pi/2;
    0 at transit and 0.5 at secondary eclipse), continuous with and within half an
    orbit of `phase`; the separation r/a = (1 - e^2)/(1 + e cos f); and
    rv = -(cos(f + omega) + e cos(omega)) / sqrt(1 - e^2), so that the planet velocity
    is Kp rv with Kp = Kp_max sin(i).
    """
    phase = np.asarray(phase, dtype=float)
    omega = np.radians(omega_deg)
    m0 = transit_mean_anomaly(e, omega_deg)
    M = 2 * np.pi * phase + m0
    E = solve_kepler(M, e)
    beta = e / (1 + np.sqrt(1 - e ** 2))
    f = E + 2 * np.arctan(beta * np.sin(E) / (1 - beta * np.cos(E)))
    rv = -(np.cos(f + omega) + e * np.cos(omega)) / np.sqrt(1 - e ** 2)
    # u = 2 pi phase + (f - M) - (f - M at transit), with the continuous equation of
    # center f - M = (f - E) + e sin(E)
    eoc = (f - E) + e * np.sin(E)
    f0 = (np.pi / 2 - omega + np.pi) % (2 * np.pi) - np.pi     # true anomaly at transit
    E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f0 / 2))
    uphase = phase + (eoc - ((f0 - E0) + e * np.sin(E0))) / (2 * np.pi)
    rscale = (1 - e ** 2) / (1 + e * np.cos(f))
    return uphase, rscale, rv


def planet_rv_eccentric(phase, kp, e, omega_deg):
    """Planet line-of-sight velocity [km/s] on an eccentric orbit.

    `phase` is (t - T0)/P with T0 the transit time, kp = Kp_max sin(i) and omega the
    star's argument of periastron [deg]: v = -kp (cos(f + omega) + e cos(omega)) /
    sqrt(1 - e^2), with the true anomaly f from Kepler's equation.
    """
    return kp * eccentric_orbit(phase, e, omega_deg)[2]


def lambert_phase(phase, inclination_deg=90.):
    """Lambertian phase curve, normalized to 1 at full phase.

    f = (sin a + (pi - a) cos a) / pi with the phase angle a given by
    cos a = -sin(i) cos(2 pi phase); for an edge-on orbit a = pi |1 - 2 phase|
    (a = 0 at secondary eclipse). For an eccentric orbit pass the orbital phase
    (eccentric_orbit); the curve is geometric only (no irradiation change with r).
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
