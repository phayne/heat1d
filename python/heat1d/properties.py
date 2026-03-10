"""Thermophysical property functions for heat1d.

Implements the material property models described in
Hayne et al. (2017), Appendix A1 (Eqs. A2-A6).
"""
import numpy as np
from astropy.constants import sigma_sb

from .config import R350 as _R350

try:
    import numba
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False


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


def heatCapacity(planet, T, model="polynomial"):
    """Calculate heat capacity of regolith (temperature-dependent).

    Parameters
    ----------
    planet : planets.Planet
        Planetary constants object
    T : float or np.ndarray
        Temperature [K]
    model : str
        Heat capacity model to use:
        - ``"polynomial"`` (default): Ledlow et al. (1992) / Hemingway et al.
          (1981) 4th-order polynomial. Valid for T > ~10 K; yields negative
          values for T < 1.3 K.
        - ``"biele2022"``: Biele et al. (2022, IJTP 43:144, Eq. 24) rational
          log-log fit to Apollo lunar data. Correctly gives Debye T^3 behavior
          at low T and c_p(0) = 0. Valid from ~5 K to ~1000 K.

    Returns
    -------
    np.ndarray
        Heat capacity cp [J kg-1 K-1]
    """
    if model == "biele2022":
        return heatCapacity_biele(T)
    c = planet.cpCoeff
    return np.polyval(c, T)


# Biele et al. (2022, IJTP 43:144) Eq. 24 coefficients for lunar regolith.
# Rational function in log-log space; p1=3 is fixed (gives Debye T^3 limit).
BIELE2022_COEFFS = {
    "p1": 3.0, "p2": -54.45, "p3": 306.8, "p4": -376.6,
    "q1": -16.81, "q2": 87.32,
}


def _heatCapacity_biele_python(T, p1, p2, p3, p4, q1, q2):
    """Biele et al. (2022) heat capacity (pure Python / numpy)."""
    T = np.asarray(T, dtype=float)
    x = np.log(np.maximum(T, 1e-10))
    numer = p1 * x ** 3 + p2 * x ** 2 + p3 * x + p4
    denom = x ** 2 + q1 * x + q2
    ln_cp = numer / denom
    return np.exp(ln_cp)


if _HAS_NUMBA:
    _heatCapacity_biele_jit = numba.njit(cache=True)(
        _heatCapacity_biele_python
    )
else:
    _heatCapacity_biele_jit = _heatCapacity_biele_python


def heatCapacity_biele(T, coeffs=None):
    """Heat capacity using the Biele et al. (2022) rational log-log model.

    Implements Eq. 24 of Biele et al. (2022, Int. J. Thermophys. 43:144):

        ln(c_p) = (p1*x^3 + p2*x^2 + p3*x + p4) / (x^2 + q1*x + q2)
        x = ln(T)

    This has no poles, correctly gives c_p -> 0 as T -> 0 with Debye T^3
    behavior, and fits mean Apollo lunar sample data to < 3% for 90-1000 K.

    Parameters
    ----------
    T : float or np.ndarray
        Temperature [K]
    coeffs : dict, optional
        Override coefficients. Default: :data:`BIELE2022_COEFFS`.

    Returns
    -------
    np.ndarray
        Heat capacity cp [J kg-1 K-1]
    """
    if coeffs is None:
        coeffs = BIELE2022_COEFFS
    return _heatCapacity_biele_python(
        T, coeffs["p1"], coeffs["p2"], coeffs["p3"],
        coeffs["p4"], coeffs["q1"], coeffs["q2"],
    )


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
