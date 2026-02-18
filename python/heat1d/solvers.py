"""Numerical solvers for the 1D heat equation.

Implements explicit, Crank-Nicolson, and fully implicit finite difference
schemes as described in Hayne et al. (2017), Appendix A2 (Eqs. A13-A20).
"""
import numpy as np

try:
    import numba
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False


def _thomas_solve_python(lower, diag, upper, rhs):
    """Pure-Python Thomas algorithm (TDMA) fallback."""
    n = len(rhs)
    cp = np.empty(n - 1)
    dp = np.empty(n)
    x = np.empty(n)

    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - lower[i - 1] * cp[i - 1]
        if i < n - 1:
            cp[i] = upper[i] / m
        dp[i] = (rhs[i] - lower[i - 1] * dp[i - 1]) / m

    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    return x


if _HAS_NUMBA:
    thomas_solve = numba.njit(cache=True)(_thomas_solve_python)
else:
    thomas_solve = _thomas_solve_python

thomas_solve.__doc__ = """Solve tridiagonal system using Thomas algorithm (TDMA).

    Solves Ax = d where A is tridiagonal with sub-diagonal `lower`,
    main diagonal `diag`, and super-diagonal `upper`.

    Parameters
    ----------
    lower : np.ndarray, length n-1
        Sub-diagonal (below main diagonal)
    diag : np.ndarray, length n
        Main diagonal
    upper : np.ndarray, length n-1
        Super-diagonal (above main diagonal)
    rhs : np.ndarray, length n
        Right-hand side vector

    Returns
    -------
    np.ndarray, length n
        Solution vector

    Notes
    -----
    O(n) time complexity. Numerically stable for diagonally dominant systems,
    which is guaranteed for the heat equation discretization since the
    diagonal element = 1 + a_i + b_i > abs(a_i) + abs(b_i).
    """


def solve_explicit(T, dt, rho, cp, k, g1, g2):
    """Explicit (forward Euler) interior temperature update.

    Implements Eq. A17 from Hayne et al. (2017).

    Parameters
    ----------
    T : np.ndarray
        Temperature array. T[0] and T[-1] must be set by BCs first.
        T[1:-1] is updated in-place.
    dt : float
        Time step [s]
    rho, cp, k : np.ndarray
        Density, heat capacity, and thermal conductivity arrays.
    g1, g2 : np.ndarray
        Geometric finite-difference coefficients (p and q in the Appendix).
    """
    alpha = g1 * k[0:-2]
    beta = g2 * k[1:-1]
    T[1:-1] = T[1:-1] + dt / (rho[1:-1] * cp[1:-1]) * (
        alpha * T[0:-2] - (alpha + beta) * T[1:-1] + beta * T[2:]
    )


def solve_implicit(T, dt, rho, cp, k, g1, g2):
    """Fully implicit (backward Euler) interior temperature update.

    Solves the tridiagonal system for the interior temperatures at the
    new time level using the Thomas algorithm.

    Parameters
    ----------
    T : np.ndarray
        Temperature array. T[0] and T[-1] must be set by BCs first.
        T[1:-1] is updated in-place.
    dt : float
        Time step [s]
    rho, cp, k : np.ndarray
        Density, heat capacity, and thermal conductivity arrays.
    g1, g2 : np.ndarray
        Geometric finite-difference coefficients.
    """
    N = len(T)
    n_interior = N - 2

    if n_interior < 1:
        return

    # Diffusion coefficients for interior nodes
    ai = dt * g1 * k[0:n_interior] / (rho[1:-1] * cp[1:-1])
    bi = dt * g2 * k[1 : n_interior + 1] / (rho[1:-1] * cp[1:-1])

    # Tridiagonal system coefficients
    sub = -ai[1:]  # sub-diagonal, length n_interior-1
    main_diag = 1.0 + ai + bi  # main diagonal, length n_interior
    sup = -bi[:-1]  # super-diagonal, length n_interior-1

    # RHS with BC contributions
    rhs = T[1:-1].copy()
    rhs[0] += ai[0] * T[0]  # surface BC contribution
    rhs[-1] += bi[-1] * T[-1]  # bottom BC contribution

    # Solve and update
    T[1:-1] = thomas_solve(sub, main_diag, sup, rhs)


def solve_crank_nicolson(T, dt, rho, cp, k, g1, g2, T0_old=None, Tn_old=None):
    """Crank-Nicolson (semi-implicit) interior temperature update.

    Second-order in time. Averages explicit and implicit contributions.

    Parameters
    ----------
    T : np.ndarray
        Temperature array. T[0] and T[-1] must be set by BCs first.
        T[1:-1] is updated in-place.
    dt : float
        Time step [s]
    rho, cp, k : np.ndarray
        Density, heat capacity, and thermal conductivity arrays.
    g1, g2 : np.ndarray
        Geometric finite-difference coefficients.
    T0_old : float, optional
        Old surface temperature (before BC update), for the explicit half.
        If None, uses current T[0].
    Tn_old : float, optional
        Old bottom temperature (before BC update), for the explicit half.
        If None, uses current T[-1].
    """
    N = len(T)
    n_interior = N - 2

    if n_interior < 1:
        return

    if T0_old is None:
        T0_old = T[0]
    if Tn_old is None:
        Tn_old = T[-1]

    # Diffusion coefficients for interior nodes
    ai = dt * g1 * k[0:n_interior] / (rho[1:-1] * cp[1:-1])
    bi = dt * g2 * k[1 : n_interior + 1] / (rho[1:-1] * cp[1:-1])

    # Half-coefficients for Crank-Nicolson
    ha = 0.5 * ai
    hb = 0.5 * bi

    # LHS tridiagonal (implicit half)
    sub = -ha[1:]
    main_diag = 1.0 + ha + hb
    sup = -hb[:-1]

    # RHS: explicit half-step using OLD temperatures
    # Build a temporary array with old boundary values for the explicit part
    T_old = T.copy()
    T_old[0] = T0_old
    T_old[-1] = Tn_old

    rhs = ha * T_old[0:-2] + (1.0 - ha - hb) * T_old[1:-1] + hb * T_old[2:]

    # Add implicit BC contributions (new boundary values)
    rhs[0] += ha[0] * T[0]
    rhs[-1] += hb[-1] * T[-1]

    # Solve and update
    T[1:-1] = thomas_solve(sub, main_diag, sup, rhs)


def computeCFL(p, config):
    """Compute the CFL-limited time step from current properties.

    Parameters
    ----------
    p : Profile
        Model profile with grid and thermophysical properties.
    config : Configurator
        Model configuration (uses F, the Fourier number).

    Returns
    -------
    float
        Maximum stable time step [s]
    """
    return np.min(config.F * p.rho[:-1] * p.cp[:-1] * p.dz ** 2 / p.k[:-1])


def getTimeStep(p, day, config):
    """Calculate the initial time step.

    For explicit solver, uses CFL stability constraint.
    For implicit/CN solvers, returns ``config.dt_init`` if set,
    otherwise ``day/24`` as a reasonable starting point for
    adaptive timestepping.  When ``adaptive_tol`` is None,
    this is the fixed step size.

    Parameters
    ----------
    p : Profile
        Model profile with grid and thermophysical properties.
    day : float
        Length of one diurnal cycle [s].
    config : Configurator
        Model configuration.

    Returns
    -------
    float
        Time step [s]
    """
    if config.solver == "explicit":
        return computeCFL(p, config)
    else:
        # For implicit/CN: use dt_init if set, otherwise day/24
        return config.dt_init if config.dt_init is not None else day / 24
