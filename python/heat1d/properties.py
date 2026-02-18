"""Thermophysical property functions for heat1d.

Implements the material property models described in
Hayne et al. (2017), Appendix A1 (Eqs. A2-A6).
"""
import numpy as np
from astropy.constants import sigma_sb

from .config import R350 as _R350


def albedoVar(A0, a, b, i):
    """Calculate solar incidence angle-dependent albedo model.

    This follows the empirical fits of Keihm (1984) and Vasavada et al. (2012).

    Parameters
    ----------
    A0 : float
        Albedo at zero solar incidence angle
    a, b : float
        Coefficients
    i : float
        Solar incidence angle [rad]
    """
    return A0 + a * (i / (np.pi / 4)) ** 3 + b * (i / (np.pi / 2)) ** 8


def T_radeq(planet, lat):
    """Calculate radiative equilibrium temperature at local noontime.

    Parameters
    ----------
    planet : planets.Planet
        Planetary constants object
    lat : float
        Latitude of point of interest [rad]

    Returns
    -------
    float
        Temperature [K]
    """
    return (
        (1 - planet.albedo)
        / (sigma_sb.value * planet.emissivity)
        * planet.S
        * np.cos(lat)
    ) ** 0.25


def T_eq(planet, lat):
    """Calculate equilibrium mean temperature for rapidly rotating bodies.

    Parameters
    ----------
    planet : planets.Planet
        Planetary constants object
    lat : float
        Latitude of point of interest [rad]

    Returns
    -------
    float
        Temperature [K]
    """
    return T_radeq(planet, lat) / np.sqrt(2)


def heatCapacity(planet, T):
    """Calculate heat capacity of regolith (temperature-dependent).

    This polynomial fit is based on data from Ledlow et al. (1992) and
    Hemingway et al. (1981), and is valid for T > ~10 K.
    The formula yields *negative* (i.e. non-physical) values for T < 1.3 K.

    Parameters
    ----------
    planet : planets.Planet
        Planetary constants object
    T : float or np.ndarray
        Temperature [K]

    Returns
    -------
    np.ndarray
        Heat capacity cp [J kg-1 K-1]
    """
    c = planet.cpCoeff
    return np.polyval(c, T)


def thermCond(kc, T, R350=_R350(2.7)):
    """Calculate temperature-dependent thermal conductivity.

    Based on Mitchell and de Pater (1994) and Vasavada et al. (2012).

    Parameters
    ----------
    kc : float or np.ndarray
        Contact (phonon) conductivity [W m-1 K-1]
    T : float or np.ndarray
        Temperature [K]
    R350 : float
        Grain-size related parameter = chi / 350^3

    Returns
    -------
    np.ndarray
        Thermal conductivity K [W m-1 K-1]
    """
    return kc * (1 + R350 * T ** 3)
