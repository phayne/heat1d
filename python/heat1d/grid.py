"""Spatial grid construction for heat1d.

Implements the non-uniform finite difference grid described in
Hayne et al. (2017), Appendix A2.2 (Eqs. A31-A33).
"""
import numpy as np


def skinDepth(P, kappa):
    """Calculate thermal skin depth.

    Parameters
    ----------
    P : float
        Period (e.g., diurnal, seasonal) [s]
    kappa : float
        Thermal diffusivity = k/(rho*cp) [m2.s-1]

    Returns
    -------
    float
        Thermal skin depth [m]
    """
    return np.sqrt(kappa * P / np.pi)


def spatialGrid(zs, m, n, b):
    """Calculate the spatial grid.

    The spatial grid is non-uniform, with layer thickness increasing
    geometrically downward: dz[i] = dz[0] * r^i where r = 1 + 1/n.

    Parameters
    ----------
    zs : float
        Thermal skin depth [m]
    m : int
        Number of layers in upper skin depth [default: 10, set in Config]
    n : int
        Layer increase with depth: dz[i] = dz[i-1]*(1+1/n) [default: 5, set in Config]
    b : int
        Number of skin depths to bottom layer [default: 20, set in Config]

    Returns
    -------
    np.ndarray
        Spatial node coordinates in meters.
    """
    dz0 = zs / m
    r = 1.0 + 1.0 / n
    zmax = zs * b
    # Compute number of layers analytically from the geometric series
    N = int(np.ceil(np.log(1 + zmax * (r - 1) / dz0) / np.log(r)))
    # Layer thicknesses: geometric progression
    dz = dz0 * r ** np.arange(N)
    # Node positions: cumulative sum with z[0] = 0
    z = np.zeros(N + 1)
    z[1:] = np.cumsum(dz)
    return z
