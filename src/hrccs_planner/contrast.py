"""Planet/star flux contrast for dayside emission.

Blackbody estimates only: the planet dayside temperature from the energy balance
of Cowan & Agol (2011), T_day = T_eff (R*/a)^1/2 (1 - A_B)^1/4 (2/3 - 5 eps/12)^1/4,
and the contrast (Rp/R*)^2 B(T_day)/B(T_eff), photon-weighted over a top-hat band.
"""
import numpy as np

H_PLANCK, C_LIGHT, K_B = 6.62607015e-34, 2.99792458e8, 1.380649e-23
RSUN, RJUP, AU = 6.957e8, 7.1492e7, 1.495978707e11
G, MSUN = 6.67430e-11, 1.98847e30
_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy >= 2 renamed trapz

# band edges [m]; IGRINS H and K, and the 2MASS-like J, L, M windows
BANDS = {'J': (1.17e-6, 1.33e-6), 'H': (1.49e-6, 1.80e-6), 'K': (1.96e-6, 2.48e-6),
         'L': (3.4e-6, 4.1e-6), 'M': (4.5e-6, 5.1e-6)}


def band_planck(teff, band):
    """Photon-weighted blackbody flux averaged over a top-hat band (arbitrary
    units); `band` is a key of BANDS or a (lo, hi) tuple in metres."""
    lo, hi = BANDS[band] if isinstance(band, str) else band
    wl = np.linspace(lo, hi, 200)
    x = H_PLANCK * C_LIGHT / (wl * K_B * teff)
    return _trapz(wl / (wl ** 5 * np.expm1(np.minimum(x, 700))), wl)


def dayside_temperature(teff, a_over_rstar, bond_albedo=0.1, redistribution=0.2):
    """Cowan & Agol (2011) dayside temperature [K]; redistribution 0 (none) to 1 (full)."""
    t0 = teff / np.sqrt(a_over_rstar)
    return t0 * (1 - bond_albedo) ** 0.25 * (2 / 3. - 5 * redistribution / 12.) ** 0.25


def dayside_contrast(teff, rstar_rsun, a_over_rstar, rp_rjup, band, bond_albedo=0.1, redistribution=0.2):
    """(planet/star dayside flux ratio in `band`, T_day [K])."""
    tday = dayside_temperature(teff, a_over_rstar, bond_albedo, redistribution)
    rp_rs = rp_rjup * RJUP / (rstar_rsun * RSUN)
    return rp_rs ** 2 * band_planck(tday, band) / band_planck(teff, band), tday


def semi_major_axis(period_d, mstar_msun):
    """Semi-major axis [m] from Kepler's third law (planet mass neglected)."""
    return (G * mstar_msun * MSUN * (period_d * 86400.) ** 2 / (4 * np.pi ** 2)) ** (1 / 3.)


def kp_max(period_d, a_m):
    """Total orbital velocity 2 pi a / P [km/s]."""
    return 2 * np.pi * a_m / (period_d * 86400.) / 1e3
