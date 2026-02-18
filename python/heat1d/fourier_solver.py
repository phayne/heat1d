"""Fourier-matrix (transmission line analogy) solver for heat1d.

Solves the 1D heat diffusion equation in the frequency domain using the
multilayer matrix method. This eliminates time-stepping and equilibration
entirely, providing ~100-1000x speedup over finite-difference solvers.

The subsurface heat conduction is solved exactly in the frequency domain
(transmission matrices computed once). The nonlinear surface radiation
term (epsilon * sigma * T^4) is handled by Newton iteration in the time
domain, using the frequency-domain conductive flux at each step.

Reference: docs/Planetary Heat Diffusion Solver Specification.pdf
"""

import numpy as np
from astropy.constants import sigma_sb

from . import orbits
from .config import R350 as _R350
from .properties import albedoVar, heatCapacity, thermCond

_SIGMA = sigma_sb.value
_OVERFLOW_GUARD = 20.0  # Re(q*d) threshold for cosh/sinh overflow


def solve_mean_temperature(F_mean, J_geo, emissivity):
    """Solve for the mean surface temperature via Newton-Raphson.

    Solves:  epsilon * sigma * T^4 = F_mean + J_geo

    Parameters
    ----------
    F_mean : float
        Mean absorbed solar flux [W/m^2].
    J_geo : float
        Geothermal heat flux [W/m^2].
    emissivity : float
        Surface infrared emissivity.

    Returns
    -------
    float
        Mean surface temperature [K].
    """
    es = emissivity * _SIGMA
    rhs = F_mean + J_geo

    # Slow-rotator approximation: T_subsolar / sqrt(2)
    # Much closer to true mean T (~215 K) than radiative eq (~285 K)
    T_subsolar = ((F_mean * np.pi + J_geo) / es) ** 0.25
    T = T_subsolar / np.sqrt(2)

    for _ in range(50):
        f = es * T ** 4 - rhs
        fp = 4 * es * T ** 3
        dT = -f / fp
        T += dT
        if abs(dT) < 1e-8:
            break

    return T


def compute_equilibrium_profile(T_mean, z, kc, rho, planet, chi, J_geo,
                                J_pump=None, k_mean_eff=None):
    """Compute the static equilibrium temperature profile T_eq(z).

    Integrates dT/dz = (J_geo - J_pump(z)) / k(T) downward from the surface
    using RK4. When J_pump is None (no rectification), reduces to the simple
    geothermal gradient dT/dz = J_geo / k(T).

    When k_mean_eff is provided, uses the exact time-averaged conductivity
    instead of computing k from T_eq. This ensures the static gradient is
    consistent with the nonlinear <k(T(t))>.

    Parameters
    ----------
    T_mean : float
        Mean surface temperature [K].
    z : np.ndarray
        Depth grid nodes [m], z[0] = 0.
    kc : np.ndarray
        Contact conductivity at each depth node [W/m/K].
    rho : np.ndarray
        Density at each depth node [kg/m^3]. (unused, kept for API consistency)
    planet : object
        Planet object (for heat capacity polynomial).
    chi : float
        Radiative conductivity parameter.
    J_geo : float
        Geothermal heat flux [W/m^2].
    J_pump : np.ndarray or None
        Rectification (thermal pumping) flux at each depth [W/m^2].
        Positive J_pump means net downward heat transport, which steepens
        the subsurface temperature gradient.
    k_mean_eff : np.ndarray or None
        Time-averaged effective conductivity at each depth [W/m/K].
        When provided, overrides k(T_eq) in the gradient computation.

    Returns
    -------
    np.ndarray
        Equilibrium temperature at each grid node [K].
    """
    R350 = _R350(chi)
    N = len(z)
    T_eq = np.empty(N)
    T_eq[0] = T_mean

    for i in range(N - 1):
        dz_i = z[i + 1] - z[i]
        # Interpolate kc at midpoint
        kc_mid = 0.5 * (kc[i] + kc[i + 1])

        # RK4 integration with interpolated J_pump
        def dTdz(T_val, kc_val, frac):
            if k_mean_eff is not None:
                k_val = (1.0 - frac) * k_mean_eff[i] + frac * k_mean_eff[i + 1]
            else:
                k_val = thermCond(kc_val, T_val, R350)
            if J_pump is None:
                net_flux = J_geo
            else:
                j_local = (1.0 - frac) * J_pump[i] + frac * J_pump[i + 1]
                net_flux = J_geo - j_local
            return net_flux / k_val

        k1 = dTdz(T_eq[i], kc[i], 0.0)
        k2 = dTdz(T_eq[i] + 0.5 * dz_i * k1, kc_mid, 0.5)
        k3 = dTdz(T_eq[i] + 0.5 * dz_i * k2, kc_mid, 0.5)
        k4 = dTdz(T_eq[i] + dz_i * k3, kc[i + 1], 1.0)

        T_eq[i + 1] = T_eq[i] + (dz_i / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    return T_eq


def compute_layer_matrices(omega, dz, k_eq, kappa_eq):
    """Build transmission matrices and compute surface impedance + depth transfer.

    For a single frequency omega, builds the 2x2 transmission matrix for each
    layer and accumulates the global product from bottom to top.

    Parameters
    ----------
    omega : float
        Angular frequency [rad/s].
    dz : np.ndarray
        Layer thicknesses [m], length N_layers - 1.
    k_eq : np.ndarray
        Thermal conductivity at each grid node [W/m/K].
    kappa_eq : np.ndarray
        Thermal diffusivity at each grid node [m^2/s].

    Returns
    -------
    Z_surf : complex
        Surface thermal impedance [K m^2 / W].
    depth_ratio : np.ndarray (complex)
        Transfer function T_hat(z_j) / T_hat_surf for each grid node.
    flux_ratio : np.ndarray (complex)
        Transfer function Q_hat(z_j) / T_hat_surf for each grid node,
        where Q is the conductive heat flux. Used for rectification.
    """
    N_layers = len(dz)
    N_z = N_layers + 1  # number of grid nodes

    # Use properties at layer midpoints (average of adjacent nodes)
    k_layer = 0.5 * (k_eq[:-1] + k_eq[1:])
    kappa_layer = 0.5 * (kappa_eq[:-1] + kappa_eq[1:])

    # Wavenumber for each layer: q = sqrt(i * omega / kappa)
    # = (1+i) * sqrt(omega / (2*kappa))
    q = np.sqrt(1j * omega / kappa_layer)

    # q*d product for each layer
    qd = q * dz

    # Accumulate the product matrix P from bottom to top.
    # P_j represents the cumulative transfer from depth z_j to the bottom.
    # At each grid node j, depth_ratio[j] = P_j[0,0] / P_total[0,0].
    #
    # To prevent overflow for high frequencies and many layers, we normalize
    # P by |P00| at each step and track the cumulative log normalization.
    # Since all outputs are ratios, the normalization factors cancel or
    # appear as well-conditioned exp(negative) terms.

    # Start from identity at the bottom (below the last layer)
    P00 = 1.0 + 0j
    P01 = 0.0 + 0j
    P10 = 0.0 + 0j
    P11 = 1.0 + 0j

    # Store normalized P[0,0] and P[1,0] at each intermediate depth
    P00_at_depth = np.empty(N_z, dtype=complex)
    P10_at_depth = np.empty(N_z, dtype=complex)
    cum_log_at_depth = np.empty(N_z, dtype=float)
    P00_at_depth[N_z - 1] = P00  # bottom boundary
    P10_at_depth[N_z - 1] = P10  # bottom boundary
    cum_log = 0.0
    cum_log_at_depth[N_z - 1] = cum_log

    # Multiply from bottom layer (index N_layers-1) to top layer (index 0)
    for j in range(N_layers - 1, -1, -1):
        # Check for overflow
        if np.real(qd[j]) > _OVERFLOW_GUARD:
            # Wave is fully attenuated in this layer.
            # Use asymptotic form: cosh(x) ~ sinh(x) ~ exp(x)/2
            # M * P becomes dominated by exp(qd) terms
            exp_qd = np.exp(qd[j])
            half_exp = 0.5 * exp_qd
            M00 = half_exp
            M01 = half_exp / (k_layer[j] * q[j])
            M10 = k_layer[j] * q[j] * half_exp
            M11 = half_exp
        else:
            cqd = np.cosh(qd[j])
            sqd = np.sinh(qd[j])
            M00 = cqd
            M01 = sqd / (k_layer[j] * q[j])
            M10 = k_layer[j] * q[j] * sqd
            M11 = cqd

        # P_new = M_j * P_old
        new_P00 = M00 * P00 + M01 * P10
        new_P01 = M00 * P01 + M01 * P11
        new_P10 = M10 * P00 + M11 * P10
        new_P11 = M10 * P01 + M11 * P11

        P00, P01, P10, P11 = new_P00, new_P01, new_P10, new_P11

        # Normalize by |P00| to prevent overflow accumulation
        norm = abs(P00)
        if norm > 0:
            P00 /= norm
            P01 /= norm
            P10 /= norm
            P11 /= norm
            cum_log += np.log(norm)

        P00_at_depth[j] = P00
        P10_at_depth[j] = P10
        cum_log_at_depth[j] = cum_log

    # Surface impedance: Z = P00 / P10 (normalization cancels in ratio)
    Z_surf = P00 / P10

    # Depth ratios: account for normalization accumulated between each depth
    # and the surface.  true_P00[j] = P00_stored[j] * exp(cum_log[j]),
    # true_P00_surf = P00_stored_surf * exp(cum_log_surf), so
    # ratio = (P00_stored[j] / P00_surf) * exp(cum_log[j] - cum_log_surf).
    log_diff = cum_log_at_depth - cum_log
    scale = np.exp(np.clip(log_diff, -500, 500))
    depth_ratio = P00_at_depth * scale / P00
    flux_ratio = P10_at_depth * scale / P00

    return Z_surf, depth_ratio, flux_ratio


def compute_layer_matrices_vectorized(omega_arr, dz, k_eq, kappa_eq):
    """Vectorized transmission matrices for all frequencies simultaneously.

    Equivalent to calling compute_layer_matrices for each frequency in omega_arr,
    but uses NumPy broadcasting to eliminate the Python loop over frequencies.

    Parameters
    ----------
    omega_arr : np.ndarray
        Angular frequencies [rad/s], shape (N_freq,).
    dz : np.ndarray
        Layer thicknesses [m], length N_layers.
    k_eq : np.ndarray
        Thermal conductivity at each grid node [W/m/K], length N_z.
    kappa_eq : np.ndarray
        Thermal diffusivity at each grid node [m^2/s], length N_z.

    Returns
    -------
    Z_surf : np.ndarray (complex), shape (N_freq,)
        Surface thermal impedance for each frequency.
    depth_ratio : np.ndarray (complex), shape (N_freq, N_z)
        Temperature transfer function T_hat(z)/T_hat_surf.
    flux_ratio : np.ndarray (complex), shape (N_freq, N_z)
        Flux transfer function Q_hat(z)/T_hat_surf.
    """
    N_layers = len(dz)
    N_z = N_layers + 1
    N_freq = len(omega_arr)

    # Properties at layer midpoints
    k_layer = 0.5 * (k_eq[:-1] + k_eq[1:])
    kappa_layer = 0.5 * (kappa_eq[:-1] + kappa_eq[1:])

    # Reshape for broadcasting: (N_freq, N_layers)
    q = np.sqrt(1j * omega_arr[:, None] / kappa_layer[None, :])
    qd = q * dz[None, :]

    # Initialize global P matrix as identity for all frequencies.
    # To prevent overflow when accumulating many layers, we normalize P
    # by |P00| at each step and track cumulative log normalization.
    P00 = np.ones(N_freq, dtype=complex)
    P01 = np.zeros(N_freq, dtype=complex)
    P10 = np.zeros(N_freq, dtype=complex)
    P11 = np.ones(N_freq, dtype=complex)

    # Store normalized results and log-normalization at each depth
    P00_at_depth = np.zeros((N_freq, N_z), dtype=complex)
    P10_at_depth = np.zeros((N_freq, N_z), dtype=complex)
    cum_log_at_depth = np.zeros((N_freq, N_z), dtype=float)
    P00_at_depth[:, -1] = P00
    P10_at_depth[:, -1] = P10
    cum_log = np.zeros(N_freq, dtype=float)
    cum_log_at_depth[:, -1] = cum_log

    # Iterate from bottom layer to top
    for j in range(N_layers - 1, -1, -1):
        # Identify which frequencies trigger the overflow guard
        real_qd = np.real(qd[:, j])
        mask_large = real_qd > _OVERFLOW_GUARD

        # Normal case: full hyperbolic (use dummy 0 where mask is True)
        qd_safe = np.where(mask_large, 0, qd[:, j])
        cqd = np.cosh(qd_safe)
        sqd = np.sinh(qd_safe)

        # Large case: asymptotic cosh(x) ~ sinh(x) ~ exp(x)/2
        half_exp = 0.5 * np.exp(np.where(mask_large, qd[:, j], 0))

        # Assemble M components
        k_q = k_layer[j] * q[:, j]
        M00 = np.where(mask_large, half_exp, cqd)
        M11 = M00
        M01 = np.where(mask_large, half_exp / k_q, sqd / k_q)
        M10 = np.where(mask_large, k_q * half_exp, k_q * sqd)

        # Matrix multiply: P_new = M_j * P_old
        new_P00 = M00 * P00 + M01 * P10
        new_P01 = M00 * P01 + M01 * P11
        new_P10 = M10 * P00 + M11 * P10
        new_P11 = M10 * P01 + M11 * P11

        P00, P01, P10, P11 = new_P00, new_P01, new_P10, new_P11

        # Normalize by |P00| to prevent overflow accumulation
        norm = np.abs(P00)
        safe = norm > 0
        P00[safe] /= norm[safe]
        P01[safe] /= norm[safe]
        P10[safe] /= norm[safe]
        P11[safe] /= norm[safe]
        cum_log[safe] += np.log(norm[safe])

        P00_at_depth[:, j] = P00
        P10_at_depth[:, j] = P10
        cum_log_at_depth[:, j] = cum_log

    # Surface impedance: normalization cancels in ratio
    Z_surf = P00 / P10

    # Depth ratios: account for normalization accumulated between each
    # depth and the surface.  log_diff <= 0 for interior points.
    log_diff = cum_log_at_depth - cum_log[:, None]
    scale = np.exp(np.clip(log_diff, -500, 500))
    depth_ratio = P00_at_depth * scale / P00[:, None]
    flux_ratio = P10_at_depth * scale / P00[:, None]

    return Z_surf, depth_ratio, flux_ratio


def compute_rectification(T_surf_hat, depth_ratios, flux_ratios, k_eq, kc, T_eq, R350, N):
    """Compute the thermal pumping (rectification) flux J_pump(z).

    The nonlinear thermal conductivity k(T) = kc*(1 + R350*T^3) produces a
    net downward heat flux from the correlation of oscillating conductivity
    and temperature gradient:

        J_pump(z) = <k'(z,t) * dT'/dz(z,t)>   (time average)

    where k' = (dk/dT)*T' is the conductivity perturbation and
    dT'/dz = -Q'(z,t)/k_eq(z) is the linearized temperature gradient
    perturbation.

    Parameters
    ----------
    T_surf_hat : np.ndarray (complex)
        One-sided FFT of converged surface temperature, length N_freq.
    depth_ratios : np.ndarray (complex), shape (N_freq, N_z)
        Temperature transfer functions T_hat(z)/T_hat_surf per frequency.
    flux_ratios : np.ndarray (complex), shape (N_freq, N_z)
        Flux transfer functions Q_hat(z)/T_hat_surf per frequency.
    k_eq : np.ndarray
        Frozen thermal conductivity at each depth [W/m/K].
    kc : np.ndarray
        Contact conductivity at each depth [W/m/K].
    T_eq : np.ndarray
        Equilibrium temperature at each depth [K].
    R350 : float
        Radiative conductivity parameter = chi / 350^3.
    N : int
        Number of time steps (full DFT length).

    Returns
    -------
    J_pump : np.ndarray
        Net downward rectification flux at each depth [W/m^2].
    """
    N_freq = len(T_surf_hat)

    # dk/dT = 3 * kc * R350 * T_eq^2  at each depth
    dk_dT = 3.0 * kc * R350 * T_eq ** 2

    # Vectorized computation over all depths and harmonics.
    # T_hat_z[n, j] = T_surf_hat[n] * depth_ratios[n, j]
    # Q_hat_z[n, j] = T_surf_hat[n] * flux_ratios[n, j]
    # Only AC harmonics (n >= 1) contribute.
    T_hat_z = T_surf_hat[1:, np.newaxis] * depth_ratios[1:]   # (N_freq-1, N_z)
    Q_hat_z = T_surf_hat[1:, np.newaxis] * flux_ratios[1:]    # (N_freq-1, N_z)

    # k_hat = dk_dT * T_hat_z
    # dTdz_hat = -Q_hat_z / k_eq
    # cross_n = Re(conj(k_hat) * dTdz_hat) = -dk_dT/k_eq * Re(conj(T_hat_z) * Q_hat_z)
    cross = -dk_dT[np.newaxis, :] / k_eq[np.newaxis, :] * np.real(
        np.conj(T_hat_z) * Q_hat_z
    )  # (N_freq-1, N_z)

    # Parseval weighting: factor 2 for all harmonics except Nyquist (last one)
    weights = np.full(N_freq - 1, 2.0)
    if N % 2 == 0:
        weights[-1] = 1.0  # Nyquist harmonic

    J_pump = np.dot(weights, cross) / N ** 2  # (N_z,)

    return J_pump


def compute_rectification_exact(T_surf_hat, depth_ratios, flux_ratios,
                                k_eq, kc, T_eq, R350, N):
    """Exact time-domain thermal pumping (replaces linear perturbation).

    Instead of linearizing dk/dT at T_eq, reconstructs the full T(t,z) and
    dTdz(t,z) in the time domain, computes exact k(T(t)) = kc*(1+R350*T^3),
    and averages k(t)*dTdz(t) to get J_pump(z).

    Also returns k_mean_eff = <k(T(t))>, the true time-averaged conductivity
    (differs from k(T_eq) because <T^3> != <T>^3).

    Parameters
    ----------
    T_surf_hat : np.ndarray (complex)
        One-sided FFT of converged surface temperature, length N_freq.
    depth_ratios : np.ndarray (complex), shape (N_freq, N_z)
        Temperature transfer functions T_hat(z)/T_hat_surf per frequency.
    flux_ratios : np.ndarray (complex), shape (N_freq, N_z)
        Flux transfer functions Q_hat(z)/T_hat_surf per frequency.
    k_eq : np.ndarray
        Frozen thermal conductivity at each depth [W/m/K].
    kc : np.ndarray
        Contact conductivity at each depth [W/m/K].
    T_eq : np.ndarray
        Equilibrium temperature at each depth [K].
    R350 : float
        Radiative conductivity parameter = chi / 350^3.
    N : int
        Number of time steps (full DFT length).

    Returns
    -------
    J_pump : np.ndarray (N_z,)
        Net downward rectification flux at each depth [W/m^2].
    k_mean_eff : np.ndarray (N_z,)
        Time-averaged effective conductivity at each depth [W/m/K].
    """
    N_freq = len(T_surf_hat)
    N_z = len(T_eq)

    # Build full depth-resolved temperature spectrum (AC + DC)
    # T_hat_z[n, j] = T_surf_hat[n] * depth_ratios[n, j]  for AC (n>=1)
    # T_hat_z[0, j] = T_eq[j] * N  (DC component)
    T_hat_z = np.zeros((N_freq, N_z), dtype=complex)
    T_hat_z[0, :] = T_eq * N
    T_hat_z[1:, :] = T_surf_hat[1:, np.newaxis] * depth_ratios[1:]

    # IFFT to get absolute T(t, z) — shape (N, N_z)
    T_t = np.fft.irfft(T_hat_z, n=N, axis=0)

    # Build temperature gradient spectrum (AC only, DC gradient = 0 for pumping)
    # dTdz_hat = -Q_hat_z / k_eq  (linearized gradient perturbation)
    Q_hat_z = np.zeros((N_freq, N_z), dtype=complex)
    Q_hat_z[1:, :] = T_surf_hat[1:, np.newaxis] * flux_ratios[1:]
    dTdz_hat = np.zeros((N_freq, N_z), dtype=complex)
    dTdz_hat[1:, :] = -Q_hat_z[1:, :] / k_eq[np.newaxis, :]

    # IFFT to get gradient fluctuation dTdz(t, z) — shape (N, N_z)
    dTdz_t = np.fft.irfft(dTdz_hat, n=N, axis=0)

    # Exact nonlinear conductivity at each (t, z)
    k_t = thermCond(kc, T_t, R350)  # shape (N, N_z)

    # Time-averaged quantities
    J_pump = np.mean(k_t * dTdz_t, axis=0)      # shape (N_z,)
    k_mean_eff = np.mean(k_t, axis=0)            # shape (N_z,)

    return J_pump, k_mean_eff


def precompute_diurnal_flux(planet, lat, nsteps, dec=0, r=None, lon=0.0):
    """Pre-compute absorbed surface flux for one diurnal cycle.

    Uses the existing insolation code (orbit geometry, angle-dependent albedo).
    For bodies with significant eccentricity (e > 0.01), computes the full
    orbital trajectory over one solar day, including time-varying distance
    and declination.

    Parameters
    ----------
    planet : object
        Planet object with albedo, S, albedoCoef, day, rAU, year,
        eccentricity, obliquity, Lp attributes.
    lat : float
        Latitude [rad].
    nsteps : int
        Number of time steps per diurnal cycle.
    dec : float
        Solar declination [rad]. Default 0 (equinox). Used only for
        circular orbit approximation (ecc <= 0.01).
    r : float or None
        Heliocentric distance [AU]. Default: planet.rAU. Used only for
        circular orbit approximation (ecc <= 0.01).
    lon : float
        Observer longitude [rad]. Default 0.

    Returns
    -------
    flux : np.ndarray
        Absorbed surface flux at each time step [W/m^2].
    dt : float
        Time step [s].
    """
    dt = planet.day / nsteps
    Sabs = planet.S * (1.0 - planet.albedo)
    t = np.arange(nsteps) * dt
    a_coef = planet.albedoCoef[0]
    b_coef = planet.albedoCoef[1]

    ecc = planet.eccentricity
    if ecc > 0.05:
        # General eccentric orbit: compute full trajectory
        P_sid = orbits.siderealPeriod(planet.day, planet.year)
        M0 = orbits.meanFromTrue(lon, ecc)
        M = M0 + orbits.TWOPI * t / planet.year
        nu = orbits.trueFromMean(M, ecc)

        # Hour angle from general formula
        h = orbits.TWOPI * t / P_sid + lon - nu

        # Time-varying heliocentric distance
        a = planet.rAU
        x = a * (1 - ecc**2)
        r_t = x / (1 + ecc * np.cos(nu))

        # Time-varying solar declination
        obliq = planet.obliquity
        Lp = planet.Lp if planet.Lp is not None else 0.0
        dec_t = np.arcsin(np.sin(obliq) * np.sin(nu + Lp))

        c = orbits.cosSolarZenith(lat, dec_t, h)
        inc = np.arccos(c)
        A_var = albedoVar(planet.albedo, a_coef, b_coef, inc)
        f = (1.0 - A_var) / (1.0 - planet.albedo)
        flux = f * Sabs * (r_t / planet.rAU) ** -2 * c
    else:
        # Circular orbit approximation (original formula)
        if r is None:
            r = planet.rAU
        h = orbits.hourAngle(t, planet.day)
        c = orbits.cosSolarZenith(lat, dec, h)
        inc = np.arccos(c)
        A_var = albedoVar(planet.albedo, a_coef, b_coef, inc)
        f = (1.0 - A_var) / (1.0 - planet.albedo)
        flux = f * Sabs * (r / planet.rAU) ** -2 * c

    return flux, dt


def solve_fourier_matrix(flux_series, dt, z, dz, kc, rho, planet, J_geo, chi=2.7,
                         max_iter=40, tol=0.5, max_outer=10, outer_tol=0.1):
    """Solve the 1D heat equation using the Fourier-matrix method.

    The subsurface conduction is solved exactly via transmission matrices.
    The nonlinear surface radiation is resolved by Newton iteration in the
    time domain. An outer loop updates the frozen material properties and
    computes the thermal pumping (rectification) correction to capture the
    solid-state greenhouse effect at depth.

    Parameters
    ----------
    flux_series : np.ndarray
        Absorbed surface flux time series [W/m^2], length N.
    dt : float
        Time step [s].
    z : np.ndarray
        Depth grid nodes [m], length N_z.
    dz : np.ndarray
        Layer thicknesses [m], length N_z - 1.
    kc : np.ndarray
        Contact conductivity at each depth node [W/m/K].
    rho : np.ndarray
        Density at each depth node [kg/m^3].
    planet : object
        Planet object with emissivity, cpCoeff attributes.
    J_geo : float
        Geothermal heat flux [W/m^2].
    chi : float
        Radiative conductivity parameter.
    max_iter : int
        Maximum Newton iterations for nonlinear radiation (inner loop).
    tol : float
        Convergence tolerance for surface temperature [K] (inner loop).
    max_outer : int
        Maximum outer iterations for property/rectification updates.
    outer_tol : float
        Convergence tolerance for mean surface temperature [K] (outer loop).

    Returns
    -------
    T_all : np.ndarray, shape (N, N_z)
        Temperature at each time step and depth [K].
    """
    N = len(flux_series)
    N_z = len(z)
    period = N * dt
    R350 = _R350(chi)
    es = planet.emissivity * _SIGMA
    F_hat = np.fft.rfft(flux_series)
    N_freq = len(F_hat)  # N//2 + 1

    _T_FLOOR = 2.0  # minimum temperature [K] to avoid division by zero

    # --- Phase 1: Mean surface temperature (initial estimate) ---
    F_mean = np.mean(flux_series)
    T_mean = solve_mean_temperature(F_mean, J_geo, planet.emissivity)

    # --- Outer loop: iterate on frozen properties and rectification ---
    J_pump = None  # no rectification on first pass
    k_mean_eff = None  # no effective conductivity on first pass
    T_mean_prev = None
    T_surf = np.full(N, T_mean)
    depth_ratios = np.zeros((N_freq, N_z), dtype=complex)

    # Frequency vector (exclude DC), precomputed once
    n_arr = np.arange(1, N_freq)
    omega_arr = 2.0 * np.pi * n_arr / period

    # Circulant index array (reused each outer iteration)
    _circ_idx = (np.arange(N)[:, None] - np.arange(N)[None, :]) % N

    for _outer in range(max_outer):
        # --- Phase 2: Equilibrium profile and frozen properties ---
        T_eq = compute_equilibrium_profile(
            T_mean, z, kc, rho, planet, chi, J_geo,
            J_pump=J_pump, k_mean_eff=k_mean_eff
        )

        k_use = thermCond(kc, T_eq, R350)
        cp_eq = heatCapacity(planet, T_eq)
        kappa_use = k_use / (rho * cp_eq)

        # --- Phase 3: Vectorized transmission matrices ---
        h_r = 4.0 * es * T_mean ** 3

        Z_surfs = np.zeros(N_freq, dtype=complex)
        depth_ratios = np.zeros((N_freq, N_z), dtype=complex)
        flux_ratios = np.zeros((N_freq, N_z), dtype=complex)

        Z_ac, dr_ac, fr_ac = compute_layer_matrices_vectorized(
            omega_arr, dz, k_use, kappa_use
        )
        Z_surfs[1:] = Z_ac
        depth_ratios[1:, :] = dr_ac
        flux_ratios[1:, :] = fr_ac

        # --- Phase 4: Initial guess for Newton ---
        # First outer iteration: linearized frequency-domain solution.
        # Subsequent iterations: reuse converged T_surf (warm start),
        # which typically needs only 1-2 Newton steps to re-converge.
        if _outer == 0:
            T_hat = np.zeros((N_freq, N_z), dtype=complex)
            T_hat[0, :] = T_eq * N
            T_hat[1:, 0] = Z_surfs[1:] * F_hat[1:] / (1.0 + h_r * Z_surfs[1:])
            T_hat[1:, :] = T_hat[1:, 0][:, None] * depth_ratios[1:, :]
            T_all = np.fft.irfft(T_hat, n=N, axis=0)
            T_surf = T_all[:, 0].copy()

        # --- Phase 5: Newton solver (circulant matrix) ---
        # Build circulant conduction matrix from admittance spectrum
        Y_diag = np.zeros(N_freq, dtype=complex)
        Y_diag[1:] = 1.0 / Z_surfs[1:]
        Y_diag = np.where(np.isfinite(Y_diag), Y_diag, 0.0)
        C_row = np.fft.irfft(Y_diag, n=N)
        C = C_row[_circ_idx]

        for _iteration in range(max_iter):
            R = flux_series + J_geo - es * T_surf ** 4 - C @ T_surf
            T_safe = np.maximum(T_surf, _T_FLOOR)
            J = np.diag(4.0 * es * T_safe ** 3)
            dT = np.linalg.solve(J + C, R)
            T_surf = np.maximum(T_surf + dT, _T_FLOOR)
            if np.max(np.abs(dT)) < tol:
                break

        # Update mean T for next outer iteration
        T_mean_new = np.mean(T_surf)

        # Compute exact time-domain rectification for next outer iteration
        T_surf_hat = np.fft.rfft(T_surf)
        k_eq_rect = thermCond(kc, T_eq, R350)
        J_pump, k_mean_eff = compute_rectification_exact(
            T_surf_hat, depth_ratios, flux_ratios,
            k_eq_rect, kc, T_eq, R350, N
        )

        # Check outer convergence on mean surface temperature
        if T_mean_prev is not None and abs(T_mean_new - T_mean_prev) < outer_tol:
            T_mean = T_mean_new
            break
        T_mean_prev = T_mean_new
        T_mean = T_mean_new

    # --- Phase 6: Reconstruct full depth profile from converged surface ---
    T_eq_final = compute_equilibrium_profile(
        T_mean, z, kc, rho, planet, chi, J_geo,
        J_pump=J_pump, k_mean_eff=k_mean_eff
    )

    T_surf_hat_final = np.fft.rfft(T_surf)
    T_hat_final = np.zeros((N_freq, N_z), dtype=complex)
    T_hat_final[0, :] = T_eq_final * N
    T_hat_final[1:, :] = T_surf_hat_final[1:, None] * depth_ratios[1:, :]

    T_all = np.fft.irfft(T_hat_final, n=N, axis=0)

    return T_all
