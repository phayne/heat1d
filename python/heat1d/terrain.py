"""Sloped-surface geometry and terrain irradiance for heat1d.

Implements insolation on a tilted planar surface (Braun & Mitchell 1983)
with self-shadowing, plus the indirect flux from the surrounding flat
terrain (thermal emission and reflected sunlight) seen by the slope at
view factor sin^2(slope/2).

Conventions
-----------
- Azimuth: radians clockwise from north viewed from above
  (0 = N, pi/2 = E, pi = S, 3*pi/2 = W), matching the JPL Horizons
  ``AZ`` convention.  ``slope_az`` is the direction the tilted surface
  faces (the downslope direction).
- Slope: dip from horizontal [rad], 0 <= slope <= pi/2.
- All fluxes returned are *absorbed* fluxes [W/m^2], consistent with
  the surface boundary condition in :mod:`heat1d.boundary`, which
  treats ``Qs`` as fully absorbed.

References
----------
Braun, J. E., & Mitchell, J. C. (1983). Solar geometry for fixed and
tracking surfaces. Solar Energy, 31(5), 439-444.

Aharonson, O., & Schorghofer, N. (2006). Subsurface ice on Mars with
rough topography. JGR, 111, E11007.
"""
import numpy as np
from astropy.constants import sigma_sb

from .properties import albedoVar

sigma = sigma_sb.value


def slope_incidence_cos(cos_z, az_sun, slope, slope_az):
    """Cosine of the solar incidence angle on a tilted surface.

    mu = cos(z)*cos(s) + sin(z)*sin(s)*cos(A_sun - az_s)

    with two shadowing clamps:

    - ``mu < 0``: sun behind the slope (self-shadowing) -> 0
    - ``cos_z <= 0``: sun below the flat horizon -> 0 (a steep slope
      can have ``mu > 0`` with the sun below the horizon; the
      surrounding flat terrain blocks it)

    Parameters
    ----------
    cos_z : float or np.ndarray
        Cosine of the solar zenith angle, **unclipped** (negative when
        the sun is below the horizon).
    az_sun : float or np.ndarray
        Solar azimuth [rad], clockwise from north.
    slope : float
        Surface slope [rad].
    slope_az : float
        Slope azimuth (downslope direction) [rad], clockwise from north.

    Returns
    -------
    mu : float or np.ndarray
        Cosine of the local incidence angle, >= 0; 0 when shadowed.
    """
    cos_z = np.asarray(cos_z, dtype=float)
    sin_z = np.sqrt(np.maximum(0.0, 1.0 - cos_z**2))
    mu = cos_z * np.cos(slope) + sin_z * np.sin(slope) * np.cos(
        np.asarray(az_sun, dtype=float) - slope_az
    )
    mu = np.where(cos_z > 0.0, np.maximum(mu, 0.0), 0.0)
    return mu if mu.ndim else float(mu)


def ground_view_factor(slope):
    """View factor of the surrounding flat ground from a tilted surface.

    F_ground = sin^2(slope/2)

    Parameters
    ----------
    slope : float
        Surface slope [rad].
    """
    return np.sin(0.5 * slope) ** 2


def sky_view_factor(slope):
    """View factor of the sky from a tilted surface.

    F_sky = cos^2(slope/2) = 1 - F_ground

    Parameters
    ----------
    slope : float
        Surface slope [rad].
    """
    return np.cos(0.5 * slope) ** 2


def direct_slope_flux(planet, cos_z, az_sun, r_au, slope, slope_az):
    """Absorbed direct solar flux on a tilted surface [W/m^2].

    Q_dir = (1 - A(i_loc)) * S * (rAU/r)^2 * mu

    where ``i_loc = arccos(mu)`` is the *local* incidence angle on the
    tilted surface, used in the angle-dependent albedo model (the
    Keihm/Vasavada A(i) is a photometric property of the regolith
    relative to its own surface normal).

    Parameters
    ----------
    planet : object
        Planet with S, rAU, albedo, albedoCoef attributes.
    cos_z : float or np.ndarray
        Unclipped cosine of the solar zenith angle.
    az_sun : float or np.ndarray
        Solar azimuth [rad], clockwise from north.
    r_au : float or np.ndarray
        Heliocentric distance [AU].
    slope, slope_az : float
        Surface slope and azimuth [rad].
    """
    mu = slope_incidence_cos(cos_z, az_sun, slope, slope_az)
    i_loc = np.arccos(np.clip(mu, 0.0, 1.0))
    a, b = planet.albedoCoef
    A = albedoVar(planet.albedo, a, b, i_loc)
    return (1.0 - A) * planet.S * (planet.rAU / np.asarray(r_au)) ** 2 * mu


def flat_scattered_radiosity(planet, cos_z, r_au):
    """Reflected-solar radiosity of flat terrain [W/m^2].

    F_scat = A_flat(theta) * S * (rAU/r)^2 * max(cos_z, 0)

    where ``theta`` is the solar zenith angle (= incidence angle on
    flat ground) and A_flat the angle-dependent albedo.

    Parameters
    ----------
    planet : object
        Planet with S, rAU, albedo, albedoCoef attributes.
    cos_z : float or np.ndarray
        Unclipped cosine of the solar zenith angle.
    r_au : float or np.ndarray
        Heliocentric distance [AU].
    """
    c = np.maximum(np.asarray(cos_z, dtype=float), 0.0)
    theta = np.arccos(np.clip(c, 0.0, 1.0))
    a, b = planet.albedoCoef
    A = albedoVar(planet.albedo, a, b, theta)
    return A * planet.S * (planet.rAU / np.asarray(r_au)) ** 2 * c


def indirect_flux_series(planet, slope, T_flat, cos_z, r_au):
    """Absorbed indirect flux on a slope from surrounding flat terrain.

    Q_ind(t) = sin^2(s/2) * [ eps^2 * sigma * T_flat(t)^4
                              + (1 - A0) * F_scat(t) ]

    The thermal term is the flat-ground radiosity ``eps*sigma*T^4``
    absorbed with IR absorptivity ``eps`` (Kirchhoff), hence ``eps^2``
    (slope and ground share ``planet.emissivity``).  The reflected-solar
    term is the flat-ground reflected radiosity absorbed with the
    normal-incidence albedo ``A0`` (the grazing-incidence boost in the
    A(i) model is a directional-beam effect and does not apply to
    hemispherically diffuse illumination).

    Approximations (documented in docs/slopes.md): infinite flat plane,
    one-way coupling (the slope does not heat the ground back), zero
    sky emission (airless body).

    Parameters
    ----------
    planet : object
        Planet with S, rAU, albedo, albedoCoef, emissivity attributes.
    slope : float
        Surface slope [rad].
    T_flat : np.ndarray
        Flat-terrain surface temperature time series [K].
    cos_z : np.ndarray
        Unclipped cosine of the solar zenith angle at the same times.
    r_au : float or np.ndarray
        Heliocentric distance [AU] at the same times.

    Returns
    -------
    Q_ind : np.ndarray
        Absorbed indirect flux [W/m^2].
    """
    hss = ground_view_factor(slope)
    eps = planet.emissivity
    Q_therm = eps * eps * sigma * np.asarray(T_flat, dtype=float) ** 4
    Q_scat = (1.0 - planet.albedo) * flat_scattered_radiosity(
        planet, cos_z, r_au
    )
    return hss * (Q_therm + Q_scat)
