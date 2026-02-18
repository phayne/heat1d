"""Boundary condition functions for heat1d.

Implements the surface and bottom boundary conditions described in
Hayne et al. (2017), Appendix A1.1 and A2.1 (Eqs. A7, A12, A21-A30).
"""
try:
    import numba
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False


def _newton_surface_python(Ts, T1, T2, dz0, kc0, R350, emissivity, sigma,
                           Qs, dtsurf, max_iter):
    """Newton solver for surface energy balance (pure Python)."""
    deltaT = Ts
    for iteration in range(max_iter):
        if abs(deltaT) <= dtsurf and iteration > 0:
            break
        x = emissivity * sigma * Ts ** 3
        kT = kc0 * (1 + R350 * Ts ** 3)
        y = 0.5 * kT / dz0

        f = x * Ts - Qs - y * (-3 * Ts + 4 * T1 - T2)
        fp = (
            4 * x
            - 3 * kc0 * R350 * Ts ** 2 * 0.5
            * (4 * T1 - 3 * Ts - T2) / dz0
            + 3 * y
        )

        deltaT = -f / fp
        Ts += deltaT
    return Ts


if _HAS_NUMBA:
    _newton_surface = numba.njit(cache=True)(_newton_surface_python)
else:
    _newton_surface = _newton_surface_python


def surfTemp(p, Qs, max_iterations=100):
    """Surface temperature calculation using Newton's root-finding method.

    Solves the surface energy balance:
        emissivity * sigma * Ts^4 - Qs - K * dT/dz|_surface = 0

    Parameters
    ----------
    p : Profile object
        Contains temperature array, grid, and thermophysical properties.
    Qs : float
        Heating rate [W.m-2] (e.g., insolation and infrared heating)
    max_iterations : int
        Maximum number of Newton iterations (default: 100).
    """
    p.T[0] = _newton_surface(
        p.T[0], p.T[1], p.T[2], p.dz[0], p.kc[0], p.surface_R350,
        p.emissivity, p.config.sigma, Qs, p.config.DTSURF, max_iterations,
    )


def botTemp(p, Qb):
    """Calculate bottom layer temperature.

    Bottom layer temperature is calculated from the interior heat
    flux and the temperature of the layer above.

    Parameters
    ----------
    p : Profile object
        Contains temperature array, grid, and thermophysical properties.
    Qb : float
        Interior heat flux [W.m-2]
    """
    p.T[-1] = p.T[-2] + (Qb / p.k[-2]) * p.dz[-1]
