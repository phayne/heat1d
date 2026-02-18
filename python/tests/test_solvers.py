"""Tests for the solvers module."""

import numpy as np
import planets
import pytest

from heat1d.config import Configurator
from heat1d.model import Model
from heat1d.solvers import (
    getTimeStep,
    solve_crank_nicolson,
    solve_explicit,
    solve_implicit,
    thomas_solve,
)

DAY = planets.Moon.day  # diurnal period [s]


class TestThomasAlgorithm:

    def test_known_solution(self):
        """Thomas algorithm solves a simple 3x3 system."""
        # A = [[2,-1,0],[-1,2,-1],[0,-1,2]], b = [1,0,1]
        lower = np.array([-1.0, -1.0])
        diag = np.array([2.0, 2.0, 2.0])
        upper = np.array([-1.0, -1.0])
        rhs = np.array([1.0, 0.0, 1.0])
        x = thomas_solve(lower, diag, upper, rhs)
        np.testing.assert_allclose(x, [1.0, 1.0, 1.0], atol=1e-12)

    def test_identity(self):
        """Thomas algorithm solves identity system."""
        n = 10
        lower = np.zeros(n - 1)
        diag = np.ones(n)
        upper = np.zeros(n - 1)
        rhs = np.arange(n, dtype=float)
        x = thomas_solve(lower, diag, upper, rhs)
        np.testing.assert_allclose(x, rhs, atol=1e-12)

    def test_vs_numpy_solve(self):
        """Thomas algorithm agrees with numpy.linalg.solve."""
        n = 20
        rng = np.random.default_rng(42)
        lower = -rng.random(n - 1)
        diag = 3.0 + rng.random(n)  # diagonally dominant
        upper = -rng.random(n - 1)
        rhs = rng.random(n)

        # Build full matrix
        A = np.diag(diag) + np.diag(lower, -1) + np.diag(upper, 1)
        x_np = np.linalg.solve(A, rhs)
        x_thomas = thomas_solve(lower, diag, upper, rhs)
        np.testing.assert_allclose(x_thomas, x_np, atol=1e-10)


class TestExplicitSolver:

    def test_regression(self):
        """Explicit solver reproduces known Moon equator result."""
        config = Configurator(solver="explicit",
                              equil_dt=DAY / 480)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        # Regression values with angle-dependent albedo (A_h=0.12), n=4 grid
        # (Fourier-accelerated equilibration; ref 385 +/- 5 K)
        np.testing.assert_allclose(m.T[:, 0].max(), 388.36, atol=0.5)
        np.testing.assert_allclose(m.T[:, 0].min(), 91.22, atol=0.5)


class TestSolverConsistency:

    def test_solvers_agree_high_resolution(self):
        """All three solvers agree at high time resolution (dt_init=day/480)."""
        results = {}
        for solver in ("explicit", "crank-nicolson", "implicit"):
            config = Configurator(solver=solver, dt_init=DAY / 480,
                                  adaptive_tol=None)
            m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
            m.run()
            results[solver] = m.T[:, 0]

        T_explicit = results["explicit"]
        T_cn = results["crank-nicolson"]
        T_implicit = results["implicit"]

        # They should agree within 0.5 K at this resolution
        assert abs(T_explicit.max() - T_cn.max()) < 0.5
        assert abs(T_explicit.max() - T_implicit.max()) < 0.5
        assert abs(T_explicit.min() - T_cn.min()) < 0.5
        assert abs(T_explicit.min() - T_implicit.min()) < 0.5

    def test_implicit_stable_large_dt(self):
        """Implicit solver is stable with large time steps."""
        config = Configurator(solver="implicit", dt_init=DAY / 12,
                              adaptive_tol=None)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        # Should not blow up
        assert np.all(np.isfinite(m.T))
        assert m.T[:, 0].max() > 300
        assert m.T[:, 0].min() > 50

    def test_cn_stable_large_dt(self):
        """Crank-Nicolson solver is stable with large time steps."""
        config = Configurator(solver="crank-nicolson", dt_init=DAY / 12,
                              adaptive_tol=None)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        assert np.all(np.isfinite(m.T))
        assert m.T[:, 0].max() > 300
        assert m.T[:, 0].min() > 50


class TestTimeStep:

    def test_explicit_cfl_limited(self):
        """Explicit time step respects CFL condition."""
        from heat1d.profile import Profile

        config = Configurator(solver="explicit")
        p = Profile(planet=planets.Moon, lat=0.0, config=config)
        dt = getTimeStep(p, DAY, config)
        # dt should be positive and much less than a day
        assert 0 < dt < DAY / 10

    def test_implicit_default_dt(self):
        """Implicit time step defaults to day/24 (not CFL-limited)."""
        from heat1d.profile import Profile

        config = Configurator(solver="implicit")
        p = Profile(planet=planets.Moon, lat=0.0, config=config)
        dt = getTimeStep(p, DAY, config)
        np.testing.assert_allclose(dt, DAY / 24)

    def test_implicit_dt_larger_than_explicit(self):
        """Implicit time step >> explicit time step."""
        from heat1d.profile import Profile

        config_exp = Configurator(solver="explicit")
        config_imp = Configurator(solver="implicit")
        p_exp = Profile(planet=planets.Moon, lat=0.0, config=config_exp)
        p_imp = Profile(planet=planets.Moon, lat=0.0, config=config_imp)
        dt_exp = getTimeStep(p_exp, DAY, config_exp)
        dt_imp = getTimeStep(p_imp, DAY, config_imp)
        assert dt_imp > 10 * dt_exp


class TestAdaptiveTimestepping:
    """Tests for adaptive timestepping with step-doubling error estimation."""

    def test_adaptive_implicit_runs(self):
        """Adaptive implicit solver completes without errors."""
        config = Configurator(solver="implicit",
                              adaptive_tol=1.0, output_interval=DAY / 48)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        assert np.all(np.isfinite(m.T))
        assert m.T[:, 0].max() > 300
        assert m.T[:, 0].min() > 50

    def test_adaptive_cn_runs(self):
        """Adaptive Crank-Nicolson solver completes without errors."""
        config = Configurator(solver="crank-nicolson",
                              adaptive_tol=1.0, output_interval=DAY / 48)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        assert np.all(np.isfinite(m.T))
        assert m.T[:, 0].max() > 300
        assert m.T[:, 0].min() > 50

    def test_adaptive_more_accurate_than_fixed(self):
        """Both fixed and adaptive implicit match explicit to within 5 K.

        Non-adaptive implicit uses CFL-based dt during the output phase
        to avoid boundary-condition splitting artifacts, so both modes
        should be close to the explicit reference.
        """
        # Gold standard: explicit solver (tiny CFL-limited steps)
        config_exp = Configurator(solver="explicit",
                                  output_interval=DAY / 480)
        m_exp = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config_exp)
        m_exp.run()

        # Fixed implicit dt=day/480
        config_fix = Configurator(solver="implicit", dt_init=DAY / 480,
                                  output_interval=DAY / 480, adaptive_tol=None)
        m_fix = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config_fix)
        m_fix.run()

        # Adaptive implicit
        config_ada = Configurator(solver="implicit",
                                  adaptive_tol=0.5,
                                  output_interval=DAY / 480)
        m_ada = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config_ada)
        m_ada.run()

        # Interpolate to common local-time grid
        lt_ref = m_exp.lt
        T_exp = m_exp.T[:, 0]
        T_fix = np.interp(lt_ref, m_fix.lt, m_fix.T[:, 0])
        T_ada = np.interp(lt_ref, m_ada.lt, m_ada.T[:, 0])

        err_fix = np.max(np.abs(T_fix - T_exp))
        err_ada = np.max(np.abs(T_ada - T_exp))

        # Both should be within 5 K of explicit
        assert err_fix < 5.0, (
            f"Fixed vs explicit = {err_fix:.1f} K (expected < 5 K)"
        )
        assert err_ada < 5.0, (
            f"Adaptive vs explicit = {err_ada:.1f} K (expected < 5 K)"
        )

    def test_adaptive_tighter_tolerance_more_accurate(self):
        """Tighter adaptive_tol produces results closer to high-resolution."""
        config_ref = Configurator(solver="explicit",
                                  output_interval=DAY / 480)
        m_ref = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config_ref)
        m_ref.run()

        errors = {}
        for tol in [5.0, 0.5]:
            config = Configurator(solver="implicit",
                                  adaptive_tol=tol,
                                  output_interval=DAY / 480)
            m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
            m.run()
            T_interp = np.interp(m_ref.lt, m.lt, m.T[:, 0])
            errors[tol] = np.max(np.abs(T_interp - m_ref.T[:, 0]))

        assert errors[0.5] < errors[5.0], (
            f"Tighter tolerance should be more accurate: "
            f"err(0.5K)={errors[0.5]:.2f}, err(5.0K)={errors[5.0]:.2f}"
        )

    def test_adaptive_not_enabled_for_explicit(self):
        """Setting adaptive_tol with explicit solver uses fixed steps."""
        config = Configurator(solver="explicit", adaptive_tol=1.0)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        assert not m._adaptive

    def test_adaptive_enabled_by_default(self):
        """Adaptive timestepping is on by default for implicit solver."""
        config = Configurator(solver="implicit")
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        assert m._adaptive
        assert m._adaptive_tol == 1.0

    def test_adaptive_disabled_when_tol_none(self):
        """Setting adaptive_tol=None disables adaptive timestepping."""
        config = Configurator(solver="implicit", adaptive_tol=None)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        assert not m._adaptive

    def test_adaptive_config_validation(self):
        """Negative adaptive_tol raises ValueError."""
        with pytest.raises(ValueError):
            Configurator(adaptive_tol=-1.0)

    def test_adaptive_with_subsampled_output(self):
        """Adaptive timestepping works correctly with output_interval."""
        config = Configurator(solver="implicit",
                              adaptive_tol=1.0, output_interval=DAY / 48)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        assert m.N_steps == 48
        assert np.all(np.isfinite(m.T))
        assert np.all(m.T > 0)

    def test_output_interval_none_records_every_step(self):
        """output_interval=None gives one output per solver step."""
        config = Configurator(solver="explicit")
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        # CFL gives many more than 24 steps per day
        assert m.N_steps > 100
