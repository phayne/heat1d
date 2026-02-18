"""Tests for the Fourier-matrix (frequency-domain) solver."""

import time

import numpy as np
import planets
import pytest

from heat1d.config import Configurator, R350
from heat1d.fourier_solver import (
    compute_equilibrium_profile,
    compute_layer_matrices,
    compute_layer_matrices_vectorized,
    compute_rectification,
    precompute_diurnal_flux,
    solve_fourier_matrix,
    solve_mean_temperature,
)
from heat1d.model import Model
from heat1d.properties import heatCapacity, thermCond


class TestSolveMeanTemperature:

    def test_basic_equilibrium(self):
        """Newton-Raphson converges to correct T for known flux."""
        from astropy.constants import sigma_sb

        sigma = sigma_sb.value
        emis = 0.95
        # Known: if T = 200 K, flux = emis * sigma * 200^4
        T_ref = 200.0
        F = emis * sigma * T_ref ** 4
        T = solve_mean_temperature(F, 0.0, emis)
        np.testing.assert_allclose(T, T_ref, atol=1e-6)

    def test_with_geothermal(self):
        """Geothermal flux raises equilibrium temperature."""
        T_no_geo = solve_mean_temperature(100.0, 0.0, 0.95)
        T_with_geo = solve_mean_temperature(100.0, 0.018, 0.95)
        assert T_with_geo > T_no_geo


class TestHomogeneousMedium:
    """Analytical check: homogeneous medium with pure cosine forcing.

    For a semi-infinite homogeneous medium with thermal diffusivity kappa
    and conductivity k, the analytical surface impedance for frequency omega is:

        Z = 1 / (k * q)  where q = (1+i) * sqrt(omega / (2*kappa))

    The surface temperature response to a cosine flux perturbation should
    have amplitude |Z| * F_amplitude / |1 + h_r * Z| and a well-defined phase.
    """

    def test_surface_amplitude_and_phase(self):
        """Single harmonic in homogeneous medium matches analytical impedance."""
        # Homogeneous medium: constant k, constant rho*cp
        k_const = 0.01  # W/m/K
        rho_const = 1500.0
        cp_const = 600.0
        kappa_const = k_const / (rho_const * cp_const)

        N = 256
        period = 2.55024e6  # lunar day [s]
        dt = period / N
        omega = 2 * np.pi / period

        # Build a uniform grid
        N_z = 50
        z_max = 1.0  # m
        z = np.linspace(0, z_max, N_z)
        dz = np.diff(z)
        kc = np.full(N_z, k_const)  # no radiative conductivity (R350=0)
        rho = np.full(N_z, rho_const)
        k_eq = np.full(N_z, k_const)
        kappa_eq = np.full(N_z, kappa_const)

        # Analytical surface impedance
        q = np.sqrt(1j * omega / kappa_const)
        Z_analytical = 1.0 / (k_const * q)

        # Numerical surface impedance from compute_layer_matrices
        Z_surf, depth_ratio, flux_ratio = compute_layer_matrices(omega, dz, k_eq, kappa_eq)

        # Z_surf should agree with analytical for a thick enough medium
        np.testing.assert_allclose(abs(Z_surf), abs(Z_analytical), rtol=0.01)
        np.testing.assert_allclose(
            np.angle(Z_surf), np.angle(Z_analytical), atol=0.05
        )


class TestDepthDecay:
    """Temperature amplitude should decay as exp(-z/delta) where delta is the
    skin depth delta = sqrt(kappa * P / pi)."""

    def test_exponential_decay(self):
        """Single-frequency amplitude decays by 1/e at one skin depth."""
        k_const = 0.01
        rho_const = 1500.0
        cp_const = 600.0
        kappa_const = k_const / (rho_const * cp_const)

        period = 2.55024e6
        omega = 2 * np.pi / period
        delta = np.sqrt(kappa_const * period / np.pi)

        # Grid extends several skin depths
        N_z = 100
        z = np.linspace(0, 5 * delta, N_z)
        dz = np.diff(z)
        k_eq = np.full(N_z, k_const)
        kappa_eq = np.full(N_z, kappa_const)

        Z_surf, depth_ratio, flux_ratio = compute_layer_matrices(omega, dz, k_eq, kappa_eq)

        # Amplitude at each depth
        amplitudes = np.abs(depth_ratio)

        # Find the grid point closest to one skin depth
        idx_delta = np.argmin(np.abs(z - delta))
        ratio_at_delta = amplitudes[idx_delta] / amplitudes[0]

        # Should be ~1/e = 0.368
        np.testing.assert_allclose(ratio_at_delta, np.exp(-1), atol=0.05)


class TestVectorizedLayerMatrices:
    """Vectorized layer matrices should match scalar version exactly."""

    def test_matches_scalar(self):
        """Vectorized function gives same results as per-frequency scalar calls."""
        k_const = 0.01
        rho_const = 1500.0
        cp_const = 600.0
        kappa_const = k_const / (rho_const * cp_const)

        period = 2.55024e6
        N = 256
        N_freq = N // 2 + 1

        N_z = 50
        z = np.linspace(0, 1.0, N_z)
        dz = np.diff(z)
        k_eq = np.full(N_z, k_const)
        kappa_eq = np.full(N_z, kappa_const)

        # Scalar: call per frequency
        omega_arr = 2.0 * np.pi * np.arange(1, N_freq) / period
        Z_scalar = np.zeros(len(omega_arr), dtype=complex)
        dr_scalar = np.zeros((len(omega_arr), N_z), dtype=complex)
        fr_scalar = np.zeros((len(omega_arr), N_z), dtype=complex)
        for i, omega in enumerate(omega_arr):
            Z_scalar[i], dr_scalar[i, :], fr_scalar[i, :] = compute_layer_matrices(
                omega, dz, k_eq, kappa_eq
            )

        # Vectorized: single call
        Z_vec, dr_vec, fr_vec = compute_layer_matrices_vectorized(
            omega_arr, dz, k_eq, kappa_eq
        )

        np.testing.assert_allclose(Z_vec, Z_scalar, rtol=1e-12)
        np.testing.assert_allclose(dr_vec, dr_scalar, rtol=1e-12)
        np.testing.assert_allclose(fr_vec, fr_scalar, rtol=1e-12)


class TestEnergyConservation:
    """The converged surface temperature must satisfy the energy balance:
    mean(eps*sigma*T^4) = mean(F) + J_geo (to within discretization error)."""

    def test_mean_energy_balance(self):
        """Mean emitted flux equals mean absorbed flux + geothermal."""
        import copy
        from astropy.constants import sigma_sb

        planet = copy.copy(planets.Moon)
        planet.H = 0.06
        config = Configurator(solver="fourier-matrix")
        lat_rad = 0.0
        nsteps = 480

        flux, dt = precompute_diurnal_flux(planet, lat_rad, nsteps)

        from heat1d.profile import Profile

        p = Profile(planet=planet, lat=lat_rad, config=config)
        T_all = solve_fourier_matrix(
            flux_series=flux, dt=dt,
            z=p.z, dz=p.dz, kc=p.kc, rho=p.rho,
            planet=planet, J_geo=planet.Qb, chi=config.chi,
        )

        T_surf = T_all[:, 0]
        es = planet.emissivity * sigma_sb.value
        mean_emitted = np.mean(es * T_surf ** 4)
        mean_absorbed = np.mean(flux) + planet.Qb

        # Energy balance should be satisfied to within ~1 W/m^2
        np.testing.assert_allclose(mean_emitted, mean_absorbed, atol=1.0)


class TestConsistencyWithTimeStepping:
    """Fourier solver should agree with time-stepping solvers.

    The Newton iteration resolves the nonlinear surface radiation exactly.
    The outer property-update loop with thermal pumping correction captures
    the solid-state greenhouse effect. Remaining differences come from
    the linearized material properties.
    """

    def test_equator_consistency(self):
        """Fourier solver agrees with implicit solver within 10 K."""
        import copy

        planet = copy.copy(planets.Moon)
        planet.H = 0.06

        config_fm = Configurator(solver="fourier-matrix", output_interval=planets.Moon.day / 480)
        m_fm = Model(planet=planet, lat=0.0, ndays=1, config=config_fm)
        m_fm.run()

        config_imp = Configurator(solver="implicit", output_interval=planets.Moon.day / 480)
        m_imp = Model(planet=planet, lat=0.0, ndays=1, config=config_imp)
        m_imp.run()

        T_fm = m_fm.T[:, 0]
        T_imp = m_imp.T[:, 0]

        assert abs(T_fm.max() - T_imp.max()) < 10, (
            f"Tmax diff: {abs(T_fm.max() - T_imp.max()):.1f} K"
        )
        assert abs(T_fm.min() - T_imp.min()) < 10, (
            f"Tmin diff: {abs(T_fm.min() - T_imp.min()):.1f} K"
        )
        assert abs(T_fm.mean() - T_imp.mean()) < 5, (
            f"Tmean diff: {abs(T_fm.mean() - T_imp.mean()):.1f} K"
        )

    def test_latitude_ordering(self):
        """Peak temperature decreases with latitude (as for time-stepping)."""
        import copy

        Tmax_by_lat = {}
        for lat_deg in (0, 30, 60):
            planet = copy.copy(planets.Moon)
            planet.H = 0.06
            lat_rad = np.deg2rad(lat_deg)
            config = Configurator(solver="fourier-matrix", output_interval=planets.Moon.day / 480)
            m = Model(planet=planet, lat=lat_rad, ndays=1, config=config)
            m.run()
            Tmax_by_lat[lat_deg] = m.T[:, 0].max()

        assert Tmax_by_lat[0] > Tmax_by_lat[30] > Tmax_by_lat[60]


class TestSpeedBenchmark:
    """Fourier solver should complete quickly despite the outer loop."""

    def test_fourier_completes_quickly(self):
        """Fourier solver with outer loop completes in well under 2 seconds."""
        import copy

        planet = copy.copy(planets.Moon)
        planet.H = 0.06

        config = Configurator(solver="fourier-matrix", output_interval=planets.Moon.day / 480)
        t0 = time.perf_counter()
        m = Model(planet=planet, lat=0.0, ndays=1, config=config)
        m.run()
        elapsed = time.perf_counter() - t0

        # The outer loop (3-5 iterations) makes this slower than before,
        # but still much faster than multi-day equilibration runs.
        assert elapsed < 2.0, f"Fourier solver took {elapsed:.3f}s (should be < 2s)"


class TestRegression:
    """Lock in expected output values to catch unintended changes."""

    def test_equator_regression(self):
        """Fourier solver equator output matches stored reference values.

        With Newton iteration, outer property-update loop, and thermal
        pumping correction, the Fourier solver closely matches
        time-stepping results (~385 K max, ~94 K min, ~211 K mean).
        """
        import copy

        planet = copy.copy(planets.Moon)
        planet.H = 0.06

        config = Configurator(solver="fourier-matrix", output_interval=planets.Moon.day / 480)
        m = Model(planet=planet, lat=0.0, ndays=1, config=config)
        m.run()

        T_surf = m.T[:, 0]
        np.testing.assert_allclose(T_surf.max(), 385.4, atol=1.0)
        np.testing.assert_allclose(T_surf.min(), 94.0, atol=1.0)
        np.testing.assert_allclose(T_surf.mean(), 210.7, atol=1.0)

    def test_output_shape(self):
        """Output array has expected shape."""
        import copy

        planet = copy.copy(planets.Moon)
        planet.H = 0.06

        config = Configurator(solver="fourier-matrix", output_interval=planets.Moon.day / 480)
        m = Model(planet=planet, lat=0.0, ndays=1, config=config)
        m.run()

        assert m.T.shape[0] == 480
        assert m.T.shape[1] == m.N_z
        assert len(m.lt) == 480
        assert np.all(np.isfinite(m.T))
        assert np.all(m.T > 0)  # no negative temperatures
