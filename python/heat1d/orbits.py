"""
This module contains functions for calculating solar
angles from orbital elements.

Includes Kepler equation solver and general hour angle
computation for eccentric orbits (e.g. Mercury 3:2 resonance).
"""

# Constants
AU = 1.49598261e11 # Astronomical Unit [m]
GM = 3.96423e-14 # G*Msun [AU**3/s**2]
TWOPI = 6.283185307

import numpy as np


# -------------------------------------------------------------------
# Kepler equation and anomaly conversions
# -------------------------------------------------------------------

def solveKepler(M, ecc, tol=1e-12, maxiter=30):
    """Solve Kepler's equation E - e*sin(E) = M via Newton-Raphson.

    Parameters
    ----------
    M : float or np.ndarray
        Mean anomaly [rad].
    ecc : float
        Orbital eccentricity (0 <= ecc < 1).
    tol : float
        Convergence tolerance [rad].
    maxiter : int
        Maximum Newton iterations.

    Returns
    -------
    E : float or np.ndarray
        Eccentric anomaly [rad].
    """
    M = np.asarray(M, dtype=float)
    E = M.copy()  # initial guess
    for _ in range(maxiter):
        dE = (E - ecc * np.sin(E) - M) / (1 - ecc * np.cos(E))
        E -= dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def trueFromMean(M, ecc):
    """Convert mean anomaly to true anomaly via Kepler's equation.

    Parameters
    ----------
    M : float or np.ndarray
        Mean anomaly [rad].
    ecc : float
        Orbital eccentricity.

    Returns
    -------
    nu : float or np.ndarray
        True anomaly [rad].
    """
    if ecc < 1e-15:
        return np.asarray(M, dtype=float)
    E = solveKepler(M, ecc)
    # tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    nu = 2.0 * np.arctan2(
        np.sqrt(1 + ecc) * np.sin(E / 2),
        np.sqrt(1 - ecc) * np.cos(E / 2),
    )
    return nu


def meanFromTrue(nu, ecc):
    """Convert true anomaly to mean anomaly (inverse Kepler).

    Parameters
    ----------
    nu : float or np.ndarray
        True anomaly [rad].
    ecc : float
        Orbital eccentricity.

    Returns
    -------
    M : float or np.ndarray
        Mean anomaly [rad].
    """
    if ecc < 1e-15:
        return np.asarray(nu, dtype=float)
    nu = np.asarray(nu, dtype=float)
    # Eccentric anomaly from true anomaly
    E = 2.0 * np.arctan2(
        np.sqrt(1 - ecc) * np.sin(nu / 2),
        np.sqrt(1 + ecc) * np.cos(nu / 2),
    )
    M = E - ecc * np.sin(E)
    return M


# -------------------------------------------------------------------
# Sidereal period and general hour angle
# -------------------------------------------------------------------

def siderealPeriod(P_solar, P_orbital):
    """Compute sidereal rotation period from solar day and orbital period.

    For prograde rotation: 1/P_sid = 1/P_solar + 1/P_orbital.

    Parameters
    ----------
    P_solar : float
        Solar day [s] (planet.day).
    P_orbital : float
        Orbital period [s] (planet.year).

    Returns
    -------
    float
        Sidereal rotation period [s].
    """
    return P_solar * P_orbital / (P_solar + P_orbital)


def generalHourAngle(t, P_sid, P_orb, ecc, lon=0.0, M0=0.0):
    """Hour angle accounting for eccentric orbit and longitude.

    h(t) = 2*pi*t/P_sid + lon - nu(t)

    where nu(t) = trueFromMean(M0 + 2*pi*t/P_orb, ecc).

    Reduces to hourAngle(t, P_solar) when ecc=0 and lon=0.

    Parameters
    ----------
    t : float or np.ndarray
        Time since t=0 [s].
    P_sid : float
        Sidereal rotation period [s].
    P_orb : float
        Orbital period [s].
    ecc : float
        Orbital eccentricity.
    lon : float
        Observer longitude [rad] (0 = subsolar meridian at perihelion).
    M0 : float
        Mean anomaly at t=0 [rad].

    Returns
    -------
    h : float or np.ndarray
        Hour angle [rad].
    """
    M = M0 + TWOPI * t / P_orb
    nu = trueFromMean(M, ecc)
    h = TWOPI * t / P_sid + lon - nu
    return h


def noonMeanAnomaly(lon, P_solar, P_orb, ecc, true_anomaly=None):
    """Mean anomaly at which local noon occurs at an observer.

    Without an ephemeris, longitude needs a reference orientation; the
    convention here is that the prime meridian (lon = 0) is subsolar at
    perihelion.  This is the IAU definition for Mercury, whose 3:2
    spin-orbit resonance puts its hot poles at 0 and 180 E and its warm
    poles at 90 and 270 E.  The subsolar east longitude at time tau
    after perihelion is nu(tau) - Omega*tau, so local noon at ``lon``
    satisfies

        2*pi*(tau - tau_mean)/P_solar = nu - M  (the equation of center),

    where tau_mean is local *mean* noon, 2*pi*tau_mean/P_solar + lon = 0
    (mod 2*pi).  When the sun's apparent motion briefly reverses (Mercury
    near perihelion) there are several solutions; the one nearest mean
    noon is used.

    Parameters
    ----------
    lon : float
        Observer east longitude [rad].
    P_solar : float
        Mean solar day [s] (planet.day).
    P_orb : float
        Orbital period [s] (planet.year).
    ecc : float
        Orbital eccentricity.
    true_anomaly : float or None
        If given, noon is placed at this true anomaly [rad] instead, and
        ``lon`` is not used.

    Returns
    -------
    float
        Mean anomaly [rad] at local noon, in [0, 2*pi).
    """
    if true_anomaly is not None:
        return float(meanFromTrue(true_anomaly, ecc)) % TWOPI

    tau_mean = (-lon % TWOPI) / TWOPI * P_solar

    def G(tau):
        M = TWOPI * tau / P_orb
        nu = trueFromMean(M, ecc)
        eoc = (nu - M + np.pi) % TWOPI - np.pi  # equation of center
        return TWOPI * (tau - tau_mean) / P_solar - eoc

    # |eoc| < pi/2 for any bound orbit of interest, so true noon lies
    # within a quarter day of mean noon.  G can be non-monotonic, so
    # find every sign change on a fine grid and refine by bisection.
    tau = tau_mean + np.linspace(-0.25, 0.25, 20001) * P_solar
    g = G(tau)
    roots = []
    for i in np.flatnonzero(np.sign(g[:-1]) != np.sign(g[1:])):
        lo, hi = tau[i], tau[i + 1]
        g_lo = g[i]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            g_mid = G(mid)
            if np.sign(g_mid) == np.sign(g_lo):
                lo, g_lo = mid, g_mid
            else:
                hi = mid
        roots.append(0.5 * (lo + hi))
    exact = tau[g == 0.0]
    roots.extend(exact)
    tau_noon = min(roots, key=lambda r: abs(r - tau_mean))
    return float(TWOPI * tau_noon / P_orb) % TWOPI


# -------------------------------------------------------------------
# Existing functions (preserved for backward compatibility)
# -------------------------------------------------------------------

def orbitParams(model):
    """Update orbital parameters from current mean anomaly.

    Computes true anomaly (from Kepler's equation), heliocentric
    distance, solar declination, and angular velocity.
    """
    a = model.planet.rAU
    ecc = model.planet.eccentricity
    obliq = model.planet.obliquity
    Lp = model.planet.Lp

    # Compute true anomaly from mean anomaly via Kepler's equation
    nu = trueFromMean(model.M, ecc)
    model.nu = nu

    # Useful parameter:
    x = a*(1 - ecc**2)

    # Distance to Sun
    model.r = x/(1 + ecc*np.cos(nu))

    # Solar declination
    model.dec = np.arcsin( np.sin(obliq)*np.sin(nu+Lp) )

    # Angular velocity (kept for reference/diagnostics)
    model.nudot = model.r**-2 * np.sqrt(GM*x)

def cosSolarZenith(lat, dec, h, clip=True):
    """Cosine of the solar zenith angle.

    Parameters
    ----------
    lat : float or np.ndarray
        Latitude [rad].
    dec : float or np.ndarray
        Solar declination [rad].
    h : float or np.ndarray
        Hour angle [rad] (0 at local noon).
    clip : bool
        When True (default), clip to zero when the sun is below the
        horizon.  When False, return the raw signed cosine (needed by
        sloped-surface geometry, which must distinguish "sun below the
        flat horizon" from "sun behind the slope").
    """
    # Cosine of solar zenith angle
    x = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(h)

    if not clip:
        return x

    # Clipping function = zero when sun below horizon:
    y = 0.5*(x + np.abs(x))

    return y


def solarAzimuth(lat, dec, h):
    """Solar azimuth angle, clockwise from north.

    Computed from the sun unit vector in the local ENU frame:

        E = -cos(dec)*sin(h)
        N =  cos(lat)*sin(dec) - sin(lat)*cos(dec)*cos(h)
        A = atan2(E, N)  in [0, 2*pi)

    The atan2 form handles all quadrants and the polar/zenith
    degeneracies without special cases (at the zenith it returns 0,
    which is harmless because the azimuth term in the slope geometry
    is multiplied by sin(z) = 0).

    Parameters
    ----------
    lat : float or np.ndarray
        Latitude [rad].
    dec : float or np.ndarray
        Solar declination [rad].
    h : float or np.ndarray
        Hour angle [rad] (0 at local noon, positive in the afternoon).

    Returns
    -------
    A : float or np.ndarray
        Solar azimuth [rad], in [0, 2*pi); 0 = N, pi/2 = E.
    """
    E = -np.cos(dec) * np.sin(h)
    N = np.cos(lat) * np.sin(dec) - np.sin(lat) * np.cos(dec) * np.cos(h)
    return np.arctan2(E, N) % TWOPI

def hourAngle(t, P):

    return (TWOPI * t/P) % TWOPI
