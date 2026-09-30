"""Tests for heat1d.diagnostics: per-layer flux probes and grid-resolution checks.

Companion to test_energy_conservation.py, which verifies the model's own
solvers conserve energy; this file verifies the new diagnostic module
that reconstructs the same conductive flux, flux divergence, and column
energy budget *after the fact*, from a completed run's saved T(t, z) --
plus the grid-vs-custom-layer resolution check these tests were written
to support (see the "thin fluffy layer" investigation referenced in
CLAUDE.md / project history).
"""

import copy
import warnings

import numpy as np
from heat1d import planets
import pytest

from heat1d.config import Configurator
from heat1d.layers import DepthLayer, compute_default_properties
from heat1d.model import Model
from heat1d.profile import Profile
from heat1d.properties import heatCapacity
from heat1d.solvers import solve_explicit

from heat1d import diagnostics as diag


# ---------------------------------------------------------------------------
# conductive_flux / flux_divergence / surface_flux_3point: exact formulas
# ---------------------------------------------------------------------------

class TestConductiveFluxExact:
    """Analytic checks against hand-built temperature/property arrays."""

    def test_conductive_flux_matches_manual_formula(self, moon, default_config):
        p = Profile(planet=moon, lat=0.0, config=default_config)
        rng = np.random.default_rng(0)
        T = 200.0 + 100.0 * rng.random(p.nlayers)
        k = diag.node_conductivity(p, T)
        F_expected = k[:-1] * np.diff(T) / p.dz
        np.testing.assert_allclose(diag.conductive_flux(p, T), F_expected)

    def test_defaults_to_profile_T(self, moon, default_config):
        p = Profile(planet=moon, lat=0.0, config=default_config)
        F_default = diag.conductive_flux(p)
        F_explicit = diag.conductive_flux(p, p.T)
        np.testing.assert_allclose(F_default, F_explicit)

    def test_flux_divergence_zero_for_uniform_temperature(self, moon, default_config):
        """No temperature gradient -> no flux anywhere -> zero divergence."""
        p = Profile(planet=moon, lat=0.0, config=default_config)
        T = np.full(p.nlayers, 250.0)
        np.testing.assert_allclose(diag.conductive_flux(p, T), 0.0)
        np.testing.assert_allclose(diag.flux_divergence(p, T), 0.0)

    def test_3point_and_2point_flux_agree_for_linear_profile_on_uniform_grid(
        self, moon, default_config,
    ):
        """A linear T(z) has a constant true gradient everywhere: on a
        UNIFORM grid, the second-order-accurate 3-point boundary
        stencil and the naive 2-point interior stencil both recover it
        exactly. (On the model's actual *non-uniform*, geometrically
        growing grid, the 3-point stencil is only exact for a linear
        profile when dz0 == dz1 -- Eq A24 implicitly assumes uniform
        spacing near the surface -- so this test forces a uniform
        near-surface grid to isolate that textbook property from the
        model's own grid non-uniformity.)"""
        p = Profile(planet=moon, lat=0.0, config=default_config)
        # Force nodes 0, 1, 2 to be uniformly spaced (the 3-point
        # stencil spans exactly these three nodes).
        dz0 = p.dz[0]
        p.z[1] = p.z[0] + dz0
        p.z[2] = p.z[1] + dz0
        p.dz[0] = dz0
        p.dz[1] = dz0
        slope = 137.0  # K/m
        T = 200.0 + slope * p.z
        F = diag.conductive_flux(p, T)
        k0 = diag.node_conductivity(p, T[:1])[0]
        np.testing.assert_allclose(F[0], k0 * slope, rtol=1e-12)
        np.testing.assert_allclose(diag.surface_flux_3point(p, T), k0 * slope,
                                   rtol=1e-12)
        assert abs(diag.bc_flux_mismatch(p, T)) < 1e-9

    def test_3point_stencil_needs_uniform_spacing_for_exactness(
        self, moon, default_config,
    ):
        """On the model's actual non-uniform grid, Eq A24's 3-point
        derivative is NOT exact even for a perfectly linear profile --
        it picks up a correction of slope*(3*dz0-dz1)/(2*dz0) instead
        of the true slope. Pinning this closed form documents a real,
        pre-existing (not custom-layer-specific) accuracy limitation of
        the surface boundary condition on a growing grid."""
        p = Profile(planet=moon, lat=0.0, config=default_config)
        slope = 137.0
        T = 200.0 + slope * p.z
        k0 = diag.node_conductivity(p, T[:1])[0]
        expected_deriv = slope * (3 * p.dz[0] - p.dz[1]) / (2 * p.dz[0])
        np.testing.assert_allclose(diag.surface_flux_3point(p, T),
                                   k0 * expected_deriv, rtol=1e-12)

    def test_3point_vs_2point_mismatch_for_quadratic_profile(self, moon):
        """For T = a + b*z + c*z^2 on a UNIFORM grid (so the 3-point
        formula is exact, but the 2-point secant is not), the mismatch
        has the closed form -k*c*dz0 (derived from the secant slope
        being the average over [0, dz0], b + c*dz0, vs. the true
        derivative b at z=0)."""
        cfg = Configurator(m=10, n=1000000)  # n huge => ~uniform first cells
        p = Profile(planet=moon, lat=0.0, config=cfg)
        b, c = 50.0, -8000.0
        T = 200.0 + b * p.z + c * p.z ** 2
        k0 = diag.node_conductivity(p, T[:1])[0]
        expected_mismatch = -k0 * c * p.dz[0]
        np.testing.assert_allclose(diag.bc_flux_mismatch(p, T),
                                   expected_mismatch, rtol=1e-3)

    def test_layer_probe_bundles_consistent_values(self, moon, default_config):
        p = Profile(planet=moon, lat=0.0, config=default_config)
        probe = diag.layer_probe(p)
        F = diag.conductive_flux(p)
        np.testing.assert_allclose(probe["F"], F)
        np.testing.assert_allclose(probe["dFdz"], diag.flux_divergence(p))
        assert probe["Qb"] == pytest.approx(F[-1])
        assert probe["bc_mismatch"] == pytest.approx(
            diag.surface_flux_3point(p) - F[0]
        )

    def test_custom_layer_node_conductivity(self, moon, default_config):
        """node_conductivity respects a per-node chi override from a
        custom layer, not just the profile-wide default."""
        layer = DepthLayer(z_top=0.0, z_bottom=0.05, rho=1100.0,
                           kc=0.002, chi=0.0, label="test")
        p = Profile(planet=moon, lat=0.0, config=default_config,
                   custom_layers=[layer])
        T = np.full(p.nlayers, 300.0)
        k = diag.node_conductivity(p, T)
        # chi=0 in the layer => no radiative term => k == kc there
        in_layer = layer.contains(p.z)
        np.testing.assert_allclose(k[in_layer], p.kc[in_layer])


# ---------------------------------------------------------------------------
# column_energy_budget: matches the solver exactly (explicit, full output)
# ---------------------------------------------------------------------------

@pytest.fixture
def explicit_full_output_model(moon):
    """Short, cheaply-equilibrated explicit run with every step saved and
    property caching disabled, so the probe can be checked bit-for-bit."""
    cfg = Configurator(solver="explicit", NYEARSEQ=1, output_interval=None)
    model = Model(planet=moon, lat=0.0, ndays=1, config=cfg)
    model.profile._prop_cache_threshold = 0.0
    model.run()
    return model


class TestColumnEnergyBudget:

    def test_requires_at_least_two_steps(self, moon, default_config):
        model = Model(planet=moon, lat=0.0, ndays=1, config=default_config)
        model.T = np.zeros((1, model.profile.nlayers))
        model.lt = np.zeros(1)
        with pytest.raises(ValueError, match="at least 2"):
            diag.column_energy_budget(model)

    def test_interior_residual_at_machine_precision(self, explicit_full_output_model):
        """Every node from index 2 onward must reproduce solve_explicit's
        own update essentially exactly -- confirms the probe's flux
        divergence formula (Eq A16) matches what the interior stencil
        actually computed, with no discrete energy-conservation defect
        anywhere except the known node-1 boundary-coupling artifact."""
        budget = diag.column_energy_budget(explicit_full_output_model)
        res = budget["residual"]
        # Six-plus orders of magnitude below node 1's boundary-artifact
        # scale (~1e2 W/m^3, see the next test) and below any physically
        # meaningful diurnal flux divergence for this run -- "machine
        # noise", not a real energy-conservation defect.
        assert np.abs(res[:, 1:]).max() < 1e-3

    def test_node_one_carries_the_known_boundary_artifact(self, explicit_full_output_model):
        """Node 1 (nearest the surface) has a real, nonzero, and LARGER
        residual than any deeper node in a homogeneous, well-resolved
        column -- surfTemp solves the new-time-level T[0] before the
        interior update runs, so node 1's explicit step mixes new T[0]
        with old T[1], T[2] (see heat1d.boundary docstring). This is an
        intrinsic baseline, not a bug; a real spatial-discretization
        problem (e.g. an unresolved custom layer) elevates MULTIPLE
        nodes above this baseline, not just node 1."""
        budget = diag.column_energy_budget(explicit_full_output_model)
        res = budget["residual"]
        node1_max = np.abs(res[:, 0]).max()
        deeper_max = np.abs(res[:, 1:]).max()
        assert node1_max > 100 * deeper_max
        assert deeper_max < 1e-3

    def test_bottom_flux_reconstructs_planet_Qb(self, explicit_full_output_model):
        model = explicit_full_output_model
        budget = diag.column_energy_budget(model)
        np.testing.assert_allclose(
            budget["Qb"], model.planet.Qb, rtol=1e-3, atol=1e-4
        )
        assert budget["Qb_expected"] == model.planet.Qb

    def test_shapes_are_consistent(self, explicit_full_output_model):
        model = explicit_full_output_model
        budget = diag.column_energy_budget(model)
        n_out = model.T.shape[0] - 1
        N = model.profile.nlayers
        assert budget["F"].shape == (n_out, N - 1)
        assert budget["dFdz"].shape == (n_out, N - 2)
        assert budget["dTdt"].shape == (n_out, N - 2)
        assert budget["residual"].shape == (n_out, N - 2)
        assert budget["bc_mismatch"].shape == (n_out,)
        assert budget["column_residual"].shape == (n_out,)
        assert len(budget["t"]) == n_out
        assert len(budget["lt"]) == n_out


# ---------------------------------------------------------------------------
# layer_resolution_check: grid-vs-layer resistance mismatch
# ---------------------------------------------------------------------------

@pytest.fixture
def europa_fluffy_layer():
    """A 2 mm, ~25x-lower-conductivity surface veneer on Europa -- the
    configuration reported to produce inaccurate results (see project
    history), used throughout this test class."""
    europa = planets.Europa
    kc_fluff = europa.ks / 25.0
    return europa, DepthLayer(z_top=0.0, z_bottom=0.002, rho=europa.rhos,
                              kc=kc_fluff, chi=2.7, label="fluffy")


class TestLayerResolutionCheck:

    def test_no_custom_layers_returns_empty(self, moon, default_config):
        p = Profile(planet=moon, lat=0.0, config=default_config)
        assert diag.layer_resolution_check(p) == []

    def test_default_grid_exactly_resolves_thin_europa_layer(self, europa_fluffy_layer):
        """Even at the default m=10, where the layer (2mm) is thinner
        than the background grid's first cell (~4.6mm), the grid now
        inserts a node exactly at the layer's bottom edge (see
        heat1d.grid.insert_custom_layer_boundaries), so the cell
        spanning the layer matches its true thickness and conductivity
        exactly. This pins the fix for the regression this module was
        built to catch (previously ratio ~2.19, a ~2x over-estimate)."""
        europa, layer = europa_fluffy_layer
        cfg = Configurator(m=10)
        p = Profile(planet=europa, lat=0.0, config=cfg, custom_layers=[layer])
        result = diag.layer_resolution_check(p)[0]
        assert result["label"] == "fluffy"
        assert result["n_nodes_inside"] == 1  # only the surface node itself
        assert len(result["cells"]) == 1
        cell = result["cells"][0]
        assert cell["cell"][2] == pytest.approx(layer.z_bottom)  # node at the edge
        assert cell["ratio"] == pytest.approx(1.0, abs=1e-9)
        assert result["max_ratio_deviation"] == pytest.approx(1.0, abs=1e-9)
        assert result["total_resistance_error"] < 1e-9

    def test_resolution_is_exact_regardless_of_grid_coarseness(self, europa_fluffy_layer):
        """total_resistance_error is ~0 at every m, not just a fine
        grid -- boundary-conforming nodes make a uniform-property
        layer's resistance exact independent of resolution elsewhere,
        unlike the raw per-cell ratio picture that mattered before the
        grid fix (a thin boundary-straddling cell used to show a large
        ratio while contributing little in absolute terms)."""
        europa, layer = europa_fluffy_layer
        for m in (10, 100, 800):
            cfg = Configurator(m=m, n=20)
            p = Profile(planet=europa, lat=0.0, config=cfg,
                       custom_layers=[layer])
            result = diag.layer_resolution_check(p)[0]
            assert result["total_resistance_error"] < 1e-9
            assert result["max_ratio_deviation"] == pytest.approx(1.0, abs=1e-9)

    def test_buried_layer_between_nodes_now_visible(self, moon, default_config):
        """A related failure mode, distinct from the reported one: a
        layer that falls strictly *between* two grid nodes -- touching
        neither -- used to be completely invisible to
        apply_custom_layers (heat1d.layers), which only overrides
        properties AT nodes. The grid fix inserts nodes at its edges
        too, so it is now applied exactly like any other layer."""
        # Build an unmodified profile first to find a gap between two
        # adjacent nodes, then place a physically significant layer
        # (same contrast as the Europa case) strictly inside that gap.
        probe_profile = Profile(planet=moon, lat=0.0, config=default_config)
        z0, z1 = probe_profile.z[5], probe_profile.z[6]
        z_top = z0 + 0.3 * (z1 - z0)
        z_bottom = z0 + 0.7 * (z1 - z0)
        kc_fluffy = moon.ks / 25.0
        layer = DepthLayer(z_top=z_top, z_bottom=z_bottom, rho=moon.rhos,
                           kc=kc_fluffy, chi=2.7, label="buried")

        p = Profile(planet=moon, lat=0.0, config=default_config,
                   custom_layers=[layer])

        # apply_custom_layers now sees it: two new nodes bound the layer,
        # and the node(s) between them carry its (not the background's)
        # conductivity.
        assert p.nlayers == probe_profile.nlayers + 2
        in_layer = layer.contains(p.z)
        assert in_layer.sum() >= 1
        np.testing.assert_allclose(p.kc[in_layer], kc_fluffy)

        result = diag.layer_resolution_check(p)[0]
        assert result["total_resistance_error"] < 1e-9
        assert result["max_ratio_deviation"] == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Profile's advisory warning (built on layer_resolution_check)
# ---------------------------------------------------------------------------

class TestLayerResolutionWarning:
    """With the grid fix (insert_custom_layer_boundaries), an ordinary
    custom layer -- of any thickness, at any depth -- is now resolved
    exactly, so this warning should essentially never fire in normal
    use. It remains as a safety net for a layer thinner than the
    merge tolerance in insert_custom_layer_boundaries (a small fraction
    of the total grid depth): such a layer gets snapped away rather
    than resolved, and should still be flagged."""

    def test_no_warning_for_ordinary_thin_layer(self, europa_fluffy_layer):
        """The reported case (2mm on Europa) at the coarsest resolution
        used to warn; it is now resolved exactly and must not warn."""
        europa, layer = europa_fluffy_layer
        cfg = Configurator(m=10)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            Profile(planet=europa, lat=0.0, config=cfg, custom_layers=[layer])

    def test_warns_for_layer_thinner_than_merge_tolerance(self, europa_fluffy_layer):
        """A layer thinner than insert_custom_layer_boundaries' snap
        tolerance (a tiny fraction of the total grid depth) collapses
        to a single point and is effectively dropped -- this remains a
        real, if extreme, gap the warning should still catch."""
        europa, layer = europa_fluffy_layer
        cfg = Configurator(m=10)
        z_max = Profile(planet=europa, lat=0.0, config=cfg).z[-1]
        eps = 1e-9 * z_max  # far below the 1e-6 * z_max merge tolerance
        vanishing = DepthLayer(z_top=0.01, z_bottom=0.01 + eps,
                               rho=layer.rho, kc=layer.kc, chi=layer.chi,
                               label="vanishing")
        with pytest.warns(UserWarning, match="not well resolved"):
            Profile(planet=europa, lat=0.0, config=cfg,
                   custom_layers=[vanishing])

    def test_no_warning_without_custom_layers(self, moon, default_config):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            Profile(planet=moon, lat=0.0, config=default_config)
