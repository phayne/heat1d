"""Tests for energy conservation in the heat1d model.

Verifies flux balance at the surface radiative boundary, the bottom
flux boundary, each interior model layer, and the full column over
a complete diurnal cycle.
"""

import copy

import numpy as np
import planets
import pytest

from heat1d.boundary import botTemp, surfTemp
from heat1d.config import Configurator
from heat1d.model import Model
from heat1d.profile import Profile
from heat1d.properties import heatCapacity, thermCond
from heat1d.solvers import solve_crank_nicolson, solve_explicit, solve_implicit


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(params=["explicit", "crank-nicolson", "implicit"])
def solver_name(request):
    return request.param


@pytest.fixture
def equilibrated_model():
    """Equator model equilibrated for 5 orbits (explicit solver)."""
    config = Configurator(solver="explicit", equil_solver="implicit",
                          NYEARSEQ=5, equil_dt=planets.Moon.day / 480)
    model = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
    model.run()
    return model


# ---------------------------------------------------------------------------
# Surface boundary energy balance
# ---------------------------------------------------------------------------

class TestSurfaceBoundaryFlux:
    """Verify the surface energy balance at each output time step."""

    def test_surface_energy_balance_each_step(self, equilibrated_model):
        """epsilon*sigma*Ts^4 = Qs + K*dT/dz at every output time step.

        The surface Newton solver should satisfy the energy balance to
        within the specified accuracy (DTSURF).  Since we don't store Qs,
        we reconstruct the balance from the stored temperature profile
        and check the residual against the conductive flux.
        """
        m = equilibrated_model
        p = m.profile
        sigma = p.config.sigma
        emis = p.emissivity
        dz0 = p.dz[0]
        kc0 = p.kc[0]
        R350 = p.R350

        residuals = []
        for i in range(len(m.T)):
            Ts = m.T[i, 0]
            T1 = m.T[i, 1]
            T2 = m.T[i, 2]

            # Emitted radiation
            rad_out = emis * sigma * Ts**4

            # Conductive flux into subsurface (second-order forward difference)
            k_surf = thermCond(kc0, Ts, R350)
            dTdz = (-3 * Ts + 4 * T1 - T2) / (2 * dz0)
            cond_flux = k_surf * dTdz

            # Energy balance: rad_out = Qs + cond_flux
            # => Qs = rad_out - cond_flux (should be >= 0 during day, = 0 at night)
            # Residual: f = rad_out - Qs - cond_flux = 0
            # Since we don't know Qs, we check that rad_out - cond_flux >= 0
            # (physically: emitted >= conducted, difference = absorbed solar)
            Qs_implied = rad_out - cond_flux
            residuals.append(Qs_implied)

        residuals = np.array(residuals)
        # Absorbed solar should be non-negative at all times
        assert np.all(residuals >= -1.0), (
            "Implied absorbed solar flux is negative "
            f"(min = {residuals.min():.2f} W/m^2)"
        )
        # At nighttime, Qs should be ~0 (radiation balances conduction)
        night_mask = (m.lt > 7.0) & (m.lt < 17.0)
        night_Qs = residuals[night_mask]
        assert np.all(np.abs(night_Qs) < 1.0), (
            "Night-time implied solar flux should be ~0, "
            f"max |Qs| = {np.abs(night_Qs).max():.2f} W/m^2"
        )

    def test_surface_newton_residual(self):
        """Newton solver residual is within DTSURF tolerance."""
        config = Configurator(solver="explicit")
        p = Profile(planet=planets.Moon, lat=0.0, config=config)

        # Apply a realistic insolation and solve
        for Qs in [0.0, 400.0, 800.0, 1200.0]:
            surfTemp(p, Qs)
            Ts = p.T[0]

            # Evaluate the energy balance residual
            k_surf = thermCond(p.kc[0], Ts, p.R350)
            dTdz = (-3 * Ts + 4 * p.T[1] - p.T[2]) / (2 * p.dz[0])
            f = p.emissivity * p.config.sigma * Ts**4 - Qs - k_surf * dTdz

            # Newton should converge to |f| within a few W/m^2
            # (DTSURF controls the temperature increment, not the flux directly)
            assert abs(f) < 5.0, (
                f"Surface residual = {f:.3f} W/m^2 at Qs = {Qs}"
            )


# ---------------------------------------------------------------------------
# Bottom boundary flux balance
# ---------------------------------------------------------------------------

class TestBottomBoundaryFlux:
    """Verify the bottom boundary condition at each output time step."""

    def test_bottom_gradient_each_step(self, equilibrated_model):
        """T[-1] = T[-2] + Qb/K[-2] * dz[-1] at every output time step."""
        m = equilibrated_model
        p = m.profile
        Qb = m.planet.Qb
        dz_last = p.dz[-1]

        for i in range(len(m.T)):
            T_bot = m.T[i, -1]
            T_above = m.T[i, -2]
            # Conductivity at second-to-last node, at that step's temperature
            k_above = thermCond(p.kc[-2], T_above, p.R350)
            expected_bot = T_above + (Qb / k_above) * dz_last

            # Allow small tolerance from the fact that k is evaluated
            # at the previous step's temperature (property update lag)
            np.testing.assert_allclose(
                T_bot, expected_bot, atol=0.5,
                err_msg=f"Bottom BC violated at step {i}: "
                f"T_bot={T_bot:.3f}, expected={expected_bot:.3f}",
            )


# ---------------------------------------------------------------------------
# Layer-by-layer conservation over a single time step
# ---------------------------------------------------------------------------

class TestSingleStepFluxBalance:
    """Verify flux balance at each interior layer for a single time step."""

    def test_explicit_step_flux_balance(self):
        """Explicit update equals conductive flux divergence at each node.

        The solver sees T[0] and T[-1] already updated by the boundary
        conditions, so we must reconstruct the solver's input array
        (new BCs + old interior) to compute the expected flux divergence.
        """
        config = Configurator(solver="explicit", NYEARSEQ=1)
        model = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)

        # Equilibrate
        while model.t < model.equiltime:
            model.advance()

        p = model.profile
        dt = model.dt

        # Record state before advance (these are the values the solver uses)
        T_old = p.T.copy()
        rho = p.rho.copy()
        cp = p.cp.copy()
        k = p.k.copy()
        g1 = p.g1
        g2 = p.g2

        # Manually apply BCs to get the solver's input
        model.updateOrbit()
        model.surfFlux()
        surfTemp(p, model.Qs)
        botTemp(p, model.planet.Qb)

        # Now p.T has new BCs but old interior: this is the solver input
        T_solver_input = p.T.copy()

        # Run the explicit solver (uses p.rho, p.cp, p.k — all still old)
        solve_explicit(p.T, dt, rho, cp, k, g1, g2)
        T_new = p.T.copy()

        # Compute expected flux divergence using the solver's actual input
        alpha = g1 * k[0:-2]
        beta = g2 * k[1:-1]
        flux_div = (
            alpha * T_solver_input[0:-2]
            - (alpha + beta) * T_solver_input[1:-1]
            + beta * T_solver_input[2:]
        )

        # Use the SAME rho, cp the solver used (pre-step values)
        dT = T_new[1:-1] - T_old[1:-1]
        energy_rate = rho[1:-1] * cp[1:-1] * dT / dt

        np.testing.assert_allclose(
            energy_rate, flux_div, rtol=1e-10,
            err_msg="Explicit step flux divergence mismatch",
        )

    def test_implicit_step_energy_balance(self):
        """Implicit step conserves total column energy (up to boundary fluxes)."""
        config = Configurator(solver="implicit", NYEARSEQ=1,
                              adaptive_tol=None)
        model = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)

        # Equilibrate
        while model.t < model.equiltime:
            model.advance()

        p = model.profile
        dt = model.dt

        # Column energy before
        E_before = np.sum(p.rho * p.cp * p.T * np.append(p.dz, p.dz[-1]))

        # Advance one step
        model.advance()

        # Column energy after
        cp_new = heatCapacity(p.planet, p.T)
        E_after = np.sum(p.rho * cp_new * p.T * np.append(p.dz, p.dz[-1]))

        # Energy change should be bounded by surface + bottom fluxes
        dE = E_after - E_before
        max_flux = model.planet.S + model.planet.Qb
        max_dE = max_flux * dt

        assert abs(dE) < max_dE, (
            f"Column energy change {dE:.2f} J/m^2 exceeds "
            f"maximum possible {max_dE:.2f} J/m^2"
        )


# ---------------------------------------------------------------------------
# Full diurnal cycle conservation
# ---------------------------------------------------------------------------

class TestDiurnalCycleConservation:
    """Verify energy conservation over a complete diurnal cycle."""

    def test_column_energy_returns_after_one_cycle(self, equilibrated_model):
        """Total stored energy at end of day matches start of day.

        For an equilibrated model, the temperature profile should be
        periodic with period = 1 day.  The total stored energy change
        over one cycle should be << the total energy throughput.
        """
        m = equilibrated_model
        p = m.profile

        # Layer thicknesses (last node gets same dz as second-to-last)
        dz_full = np.append(p.dz, p.dz[-1])

        T_first = m.T[0, :]
        T_last = m.T[-1, :]

        # Stored energy change per unit area
        cp_first = heatCapacity(p.planet, T_first)
        stored_change = np.sum(p.rho * cp_first * (T_last - T_first) * dz_full)

        # Reference: total absorbed solar energy over one day
        # S * (1-A) * 1/pi * day  (average over hemisphere)
        ref_energy = m.planet.S * (1 - m.planet.albedo) / np.pi * m.planet.day

        relative_error = abs(stored_change) / ref_energy
        assert relative_error < 0.01, (
            f"Column energy imbalance = {relative_error:.4f} "
            f"({stored_change:.2f} J/m^2 over one day)"
        )

    def test_layer_temperatures_periodic(self, equilibrated_model):
        """Each layer's temperature is periodic after equilibration.

        The surface and near-surface layers should return to nearly the
        same temperature after one full diurnal cycle.
        """
        m = equilibrated_model
        T_first = m.T[0, :]
        T_last = m.T[-1, :]

        # Surface and near-surface (first 10 layers) should be within 1 K
        n_check = min(10, len(T_first))
        np.testing.assert_allclose(
            T_first[:n_check], T_last[:n_check], atol=1.0,
            err_msg="Near-surface layers are not periodic after equilibration",
        )

        # Deep layers should be even more stable (within 0.01 K)
        if len(T_first) > 20:
            np.testing.assert_allclose(
                T_first[20:], T_last[20:], atol=0.01,
                err_msg="Deep layers are not stable after equilibration",
            )

    @pytest.mark.parametrize("solver", ["explicit", "crank-nicolson", "implicit"])
    def test_energy_conservation_all_solvers(self, solver):
        """Energy conservation holds for all three solver schemes."""
        config = Configurator(solver=solver, NYEARSEQ=3, adaptive_tol=None)
        model = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        model.run()

        p = model.profile
        dz_full = np.append(p.dz, p.dz[-1])

        T_first = model.T[0, :]
        T_last = model.T[-1, :]
        cp_first = heatCapacity(p.planet, T_first)
        stored_change = np.sum(p.rho * cp_first * (T_last - T_first) * dz_full)

        ref_energy = model.planet.S * (1 - model.planet.albedo) / np.pi * model.planet.day
        relative_error = abs(stored_change) / ref_energy

        assert relative_error < 0.01, (
            f"[{solver}] Column energy imbalance = {relative_error:.4f}"
        )

    def test_nighttime_monotonic_cooling(self, equilibrated_model):
        """Surface temperature decreases monotonically during the night.

        This is a physics consistency check: with no solar input and
        positive outgoing radiation, the surface must cool continuously.
        """
        m = equilibrated_model
        T_surf = m.T[:, 0]
        lt = m.lt

        # Deep night: 8 to 16 hours past noon (well away from sunrise/set)
        night = (lt > 8.0) & (lt < 16.0)
        T_night = T_surf[night]

        # Temperature should be monotonically decreasing
        dT = np.diff(T_night)
        assert np.all(dT <= 0), (
            f"Surface temperature is not monotonically decreasing at night "
            f"(max increase = {dT.max():.4f} K)"
        )

    def test_subsurface_temperature_exceeds_surface_mean(self, equilibrated_model):
        """Mean subsurface T > mean surface T (Jensen's inequality).

        The nonlinear T^4 radiation law means the surface radiates more
        efficiently during the day than it absorbs at night, depressing
        the surface mean below the subsurface mean.
        """
        m = equilibrated_model
        mean_surf = m.T[:, 0].mean()

        # Check a few subsurface layers
        for j in [5, 10, 15]:
            if j < m.T.shape[1]:
                mean_sub = m.T[:, j].mean()
                assert mean_sub > mean_surf, (
                    f"Layer {j}: mean subsurface T ({mean_sub:.2f} K) "
                    f"should exceed mean surface T ({mean_surf:.2f} K)"
                )
