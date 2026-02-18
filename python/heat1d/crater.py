"""Bowl-shaped crater PSR physics (Ingersoll & Svitek 1992).

Implements the analytical framework for computing surface temperatures
in permanently shadowed regions (PSRs) of bowl-shaped (spherical cap)
craters.  The crater floor receives no direct sunlight but is heated
by scattered sunlight and thermal IR emitted by the sunlit crater walls.

Reference
---------
Ingersoll, A. P., Svitek, T., & Murray, B. C. (1992).
Stability of polar frosts in spherical bowl-shaped craters on the Moon,
Mercury, and Mars. *Icarus*, 100, 40-47.
"""

import numpy as np


def crater_f(d_D):
    """Surface area ratio *f* from crater depth/diameter ratio.

    Parameters
    ----------
    d_D : float
        Crater depth-to-diameter ratio (typically 0.1-0.2).

    Returns
    -------
    float
        Surface area ratio f = (1 - cos beta) / 2.
    """
    return 4.0 * d_D**2 / (1.0 + 4.0 * d_D**2)


def crater_beta(f):
    """Crater half-angle beta from surface area ratio.

    Parameters
    ----------
    f : float
        Surface area ratio.

    Returns
    -------
    float
        Half-angle beta [rad].
    """
    return np.arccos(1.0 - 2.0 * f)


def psr_viable(lat_rad, obliquity_rad, beta):
    """Check whether a PSR can exist at the given latitude.

    A permanently shadowed region exists when the maximum solar
    elevation angle is less than the crater half-angle beta.

    Parameters
    ----------
    lat_rad : float
        Latitude [rad].
    obliquity_rad : float
        Obliquity of the spin axis [rad].
    beta : float
        Crater half-angle [rad].

    Returns
    -------
    viable : bool
        True if PSR can exist.
    e0_max : float
        Maximum solar elevation angle [rad].
    """
    e0_max = np.pi / 2.0 - abs(lat_rad) + obliquity_rad
    return e0_max < beta, e0_max


def effective_emissivity(emissivity, f):
    """Effective emissivity for a cavity (Ingersoll Eq. 7).

    epsilon_eff = epsilon / [1 - (1 - epsilon) * f]

    Parameters
    ----------
    emissivity : float
        Flat-surface IR emissivity.
    f : float
        Crater surface area ratio.

    Returns
    -------
    float
        Effective cavity emissivity (>= emissivity).
    """
    return emissivity / (1.0 - (1.0 - emissivity) * f)


def psr_flux(F0, sin_e0, f, albedo, emissivity):
    """Absorbed flux at the PSR crater floor (Ingersoll Eq. 8).

    Q_psr = F0 * sin(e0) * f * (1-A) / (1-A*f) * [eps + A*(1-f)]

    This is the effective Q_s to use with the *original* emissivity
    in the surface energy balance.  The cavity IR trapping is handled
    separately via ``effective_emissivity()``.

    Parameters
    ----------
    F0 : float or array
        Top-of-atmosphere solar flux [W/m^2] (= S / R_AU^2).
    sin_e0 : float or array
        Sine of the solar elevation angle (clamped >= 0).
    f : float
        Crater surface area ratio.
    albedo : float
        Bond albedo (A_0).
    emissivity : float
        Flat-surface IR emissivity.

    Returns
    -------
    float or array
        Absorbed flux Q_psr [W/m^2].
    """
    return (F0 * sin_e0
            * f * (1.0 - albedo) / (1.0 - albedo * f)
            * (emissivity + albedo * (1.0 - f)))
