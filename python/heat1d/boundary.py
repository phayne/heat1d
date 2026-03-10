"""Boundary condition functions for heat1d.

Implements the surface and bottom boundary conditions described in
Hayne et al. (2017), Appendix A1.1 and A2.1 (Eqs. A7, A12, A21-A30).

Stability Notes
---------------
The surface energy balance involves a nonlinear Stefan-Boltzmann T^4
radiation term. Williams & Curry (1977, IJNME 11:1605) showed that
standard implicit finite-difference methods can produce spurious surface
temperature oscillations when the T^4 boundary condition is linearized
and embedded directly into the tridiagonal system.

This module avoids that instability through **operator splitting**: the
surface temperature is solved to full nonlinear convergence via Newton's
method *before* the interior solver runs. The interior solver then sees
fixed Dirichlet boundary values and remains unconditionally stable. This
contrasts with linearization approaches (e.g., Schorghofer & Khatiwala
2024, PSJ 5:120) that approximate T^4 around a reference temperature
and fold it into the tridiagonal matrix, which can go unstable when the
reference temperature is a poor estimate (e.g., at sunrise).

The trade-off is that operator splitting is first-order accurate in the
boundary-condition coupling (vs. second-order for a fully coupled
Crank-Nicolson scheme), but it avoids the conditional instability
entirely.
"""
import math

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


def _volterra_predictor_python(T0, T1, T2, Qs_prev, Qs_new, dt,
                               kc0, R350, rho0, cp0, emissivity, sigma, dz0):
    """Volterra integral predictor for surface temperature (pure Python).

    Provides a second-order accurate initial guess for the Newton
    iteration, following Schorghofer & Khatiwala (2024, PSJ 5:120, Eq. 62).

    Parameters
    ----------
    T0 : float
        Surface temperature at the previous time step [K].
    T1, T2 : float
        Temperatures at nodes 1 and 2 (for conductive flux estimate).
    Qs_prev, Qs_new : float
        Absorbed solar flux at previous and current time steps [W/m^2].
    dt : float
        Time step [s].
    kc0 : float
        Contact conductivity at the surface [W/m/K].
    R350 : float
        chi/350^3 radiative conductivity parameter.
    rho0, cp0 : float
        Density and heat capacity at the surface.
    emissivity, sigma : float
        Emissivity and Stefan-Boltzmann constant.
    dz0 : float
        First grid spacing [m].

    Returns
    -------
    float
        Predicted surface temperature [K].
    """
    k0 = kc0 * (1.0 + R350 * T0 ** 3)
    # Conductive heat flux at surface (3-point stencil, same as Newton solver)
    H0 = k0 * (-3.0 * T0 + 4.0 * T1 - T2) / (2.0 * dz0)
    # Thermal inertia
    Gamma = math.sqrt(rho0 * cp0 * k0)
    # Weighted flux average (linear interpolation of Q over [0, dt])
    Q_avg = (2.0 * Qs_new + Qs_prev) / 3.0
    # Radiation at previous time
    rad = emissivity * sigma * T0 ** 4
    # Denominator (diffusive + radiative damping)
    denom = math.sqrt(math.pi / (4.0 * dt)) * Gamma + (8.0 / 3.0) * emissivity * sigma * T0 ** 3
    if denom < 1e-30:
        return T0
    delta_T = (-H0 - rad + Q_avg) / denom
    T_pred = T0 + delta_T
    return max(T_pred, 2.0)


if _HAS_NUMBA:
    _volterra_predictor = numba.njit(cache=True)(_volterra_predictor_python)
else:
    _volterra_predictor = _volterra_predictor_python


def surfTemp(p, Qs, max_iterations=100, Qs_prev=None, dt=None):
    """Surface temperature calculation using Newton's root-finding method.

    Solves the surface energy balance:
        emissivity * sigma * Ts^4 - Qs - K * dT/dz|_surface = 0

    When ``p.config.use_volterra_predictor`` is True and both *Qs_prev*
    and *dt* are provided, the Volterra integral predictor (Schorghofer
    & Khatiwala 2024, Eq. 62) is used to generate a better initial guess
    for the Newton iteration.

    Parameters
    ----------
    p : Profile object
        Contains temperature array, grid, and thermophysical properties.
    Qs : float
        Heating rate [W.m-2] (e.g., insolation and infrared heating)
    max_iterations : int
        Maximum number of Newton iterations (default: 100).
    Qs_prev : float, optional
        Absorbed flux at the previous time step [W/m^2].
    dt : float, optional
        Current time step [s].
    """
    if (p.config.use_volterra_predictor
            and Qs_prev is not None and dt is not None and dt > 0):
        p.T[0] = _volterra_predictor(
            p.T[0], p.T[1], p.T[2], Qs_prev, Qs, dt,
            p.kc[0], p.surface_R350, p.rho[0], p.cp[0],
            p.emissivity, p.config.sigma, p.dz[0],
        )
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
