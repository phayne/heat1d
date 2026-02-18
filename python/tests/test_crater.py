"""Tests for bowl-shaped crater PSR physics (Ingersoll & Svitek 1992)."""

import numpy as np
import pytest

from heat1d.crater import (
    crater_beta,
    crater_f,
    effective_emissivity,
    psr_flux,
    psr_viable,
)


# ---- crater_f ----

class TestCraterF:
    def test_d_D_half(self):
        """d/D = 0.5 should give f = 0.5."""
        assert crater_f(0.5) == pytest.approx(0.5)

    def test_d_D_zero(self):
        """d/D = 0 (flat) should give f = 0."""
        assert crater_f(0.0) == pytest.approx(0.0)

    def test_d_D_02(self):
        """d/D = 0.2 should give f = 4*0.04/(1+0.16) = 0.16/1.16."""
        expected = 0.16 / 1.16
        assert crater_f(0.2) == pytest.approx(expected)

    def test_large_d_D(self):
        """As d/D → ∞, f → 1."""
        assert crater_f(100.0) == pytest.approx(1.0, abs=1e-4)


# ---- crater_beta ----

class TestCraterBeta:
    def test_f_zero(self):
        """f = 0 → β = 0."""
        assert crater_beta(0.0) == pytest.approx(0.0)

    def test_f_half(self):
        """f = 0.5 → β = π/2."""
        assert crater_beta(0.5) == pytest.approx(np.pi / 2)

    def test_f_one(self):
        """f = 1 → β = π."""
        assert crater_beta(1.0) == pytest.approx(np.pi)

    def test_round_trip(self):
        """crater_beta(crater_f(d_D)) should give the expected β."""
        d_D = 0.2
        f = crater_f(d_D)
        beta = crater_beta(f)
        # Verify f = (1 - cos β) / 2
        assert f == pytest.approx((1.0 - np.cos(beta)) / 2.0)


# ---- psr_viable ----

class TestPSRViable:
    def test_polar(self):
        """PSR at 89°N on the Moon (obliquity ~1.54°) should be viable."""
        lat = np.deg2rad(89.0)
        obliquity = np.deg2rad(1.54)
        f = crater_f(0.2)
        beta = crater_beta(f)
        viable, e0_max = psr_viable(lat, obliquity, beta)
        assert viable
        assert e0_max < beta

    def test_equator(self):
        """PSR at equator should NOT be viable for any reasonable d/D."""
        lat = 0.0
        obliquity = np.deg2rad(1.54)
        f = crater_f(0.2)
        beta = crater_beta(f)
        viable, e0_max = psr_viable(lat, obliquity, beta)
        assert not viable
        assert e0_max > beta

    def test_e0_max_value(self):
        """Verify e0_max = π/2 - |lat| + obliquity."""
        lat = np.deg2rad(80.0)
        obliquity = np.deg2rad(1.54)
        f = crater_f(0.2)
        beta = crater_beta(f)
        _, e0_max = psr_viable(lat, obliquity, beta)
        expected = np.pi / 2 - np.deg2rad(80.0) + np.deg2rad(1.54)
        assert e0_max == pytest.approx(expected)

    def test_southern_hemisphere(self):
        """Negative latitude should also work (uses abs)."""
        lat = np.deg2rad(-89.0)
        obliquity = np.deg2rad(1.54)
        f = crater_f(0.2)
        beta = crater_beta(f)
        viable, _ = psr_viable(lat, obliquity, beta)
        assert viable


# ---- effective_emissivity ----

class TestEffectiveEmissivity:
    def test_f_zero(self):
        """f = 0 → ε_eff = ε (flat surface, no cavity effect)."""
        assert effective_emissivity(0.95, 0.0) == pytest.approx(0.95)

    def test_increases_with_f(self):
        """Effective emissivity should increase with f (cavity trapping)."""
        eps = 0.95
        e1 = effective_emissivity(eps, 0.1)
        e2 = effective_emissivity(eps, 0.3)
        assert e1 > eps
        assert e2 > e1

    def test_perfect_emitter(self):
        """ε = 1.0 → ε_eff = 1.0 regardless of f."""
        assert effective_emissivity(1.0, 0.5) == pytest.approx(1.0)

    def test_known_value(self):
        """Verify against hand calculation: ε=0.95, f=0.138."""
        f = crater_f(0.2)
        eps_eff = effective_emissivity(0.95, f)
        expected = 0.95 / (1.0 - 0.05 * f)
        assert eps_eff == pytest.approx(expected)


# ---- psr_flux ----

class TestPSRFlux:
    def test_nighttime(self):
        """sin(e₀) = 0 → Q_psr = 0."""
        assert psr_flux(1361.0, 0.0, 0.138, 0.12, 0.95) == pytest.approx(0.0)

    def test_zero_albedo(self):
        """A = 0: Q_psr = F₀ sin(e₀) · f · ε."""
        F0 = 1361.0
        sin_e0 = 0.5
        f = 0.2
        A = 0.0
        eps = 0.95
        expected = F0 * sin_e0 * f * eps
        assert psr_flux(F0, sin_e0, f, A, eps) == pytest.approx(expected)

    def test_positive(self):
        """PSR flux should be positive when sun is up."""
        q = psr_flux(1361.0, 0.1, 0.138, 0.12, 0.95)
        assert q > 0

    def test_less_than_direct(self):
        """PSR flux should be much less than direct illumination."""
        F0 = 1361.0
        sin_e0 = 1.0  # overhead sun
        f = crater_f(0.2)
        A = 0.12
        eps = 0.95
        q_psr = psr_flux(F0, sin_e0, f, A, eps)
        q_direct = F0 * (1 - A) * sin_e0
        assert q_psr < q_direct / 3  # PSR flux << direct

    def test_array_input(self):
        """Should work with numpy arrays."""
        F0 = np.array([1361.0, 1361.0, 0.0])
        sin_e0 = np.array([0.5, 0.0, 0.3])
        q = psr_flux(F0, sin_e0, 0.138, 0.12, 0.95)
        assert q[1] == pytest.approx(0.0)
        assert q[2] == pytest.approx(0.0)
        assert q[0] > 0


# ---- Integration test with Model ----

class TestPSRModel:
    def test_psr_gives_cold_temperature(self):
        """PSR crater at high latitude should produce much colder temperatures."""
        import planets
        from heat1d.config import Configurator
        from heat1d.model import Model

        config = Configurator(solver="implicit", NYEARSEQ=1, m=10, n=4, b=15)
        config.output_interval = planets.Moon.day / 48

        # Flat surface at 85°N
        model_flat = Model(
            planet=planets.Moon,
            lat=np.deg2rad(85.0),
            ndays=1,
            config=config,
        )
        model_flat.run()

        # PSR crater at 85°N
        model_psr = Model(
            planet=planets.Moon,
            lat=np.deg2rad(85.0),
            ndays=1,
            config=config,
            psr_d_D=0.2,
        )
        model_psr.run()

        T_max_flat = model_flat.T[:, 0].max()
        T_max_psr = model_psr.T[:, 0].max()

        # PSR should be significantly colder
        assert T_max_psr < T_max_flat - 20

    def test_psr_warning_at_equator(self):
        """PSR at equator should produce a warning."""
        import warnings
        import planets
        from heat1d.config import Configurator
        from heat1d.model import Model

        config = Configurator(solver="implicit", NYEARSEQ=1, m=10, n=4, b=10)
        config.output_interval = planets.Moon.day / 24

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            model = Model(
                planet=planets.Moon,
                lat=0.0,
                ndays=1,
                config=config,
                psr_d_D=0.2,
            )
            assert len(w) == 1
            assert "PSR not viable" in str(w[0].message)
