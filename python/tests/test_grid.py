"""Tests for the grid module."""

import numpy as np
from heat1d import planets
import pytest

from heat1d.grid import (
    insert_custom_layer_boundaries, skinDepth, spatialGrid,
)
from heat1d.layers import DepthLayer


class TestSkinDepth:

    def test_known_value(self):
        """Skin depth for known values."""
        P = 2.55e6  # Moon synodic day [s]
        kappa = 7.4e-4 / (1100.0 * 600.0)  # ks / (rhos * cp0)
        zs = skinDepth(P, kappa)
        # Should be a few cm
        assert 0.01 < zs < 0.1

    def test_proportional_to_sqrt_period(self):
        """Skin depth scales as sqrt(P)."""
        kappa = 1e-6
        zs1 = skinDepth(1.0, kappa)
        zs4 = skinDepth(4.0, kappa)
        np.testing.assert_allclose(zs4 / zs1, 2.0, rtol=1e-10)

    def test_proportional_to_sqrt_kappa(self):
        """Skin depth scales as sqrt(kappa)."""
        P = 1.0
        zs1 = skinDepth(P, 1.0)
        zs4 = skinDepth(P, 4.0)
        np.testing.assert_allclose(zs4 / zs1, 2.0, rtol=1e-10)


class TestSpatialGrid:

    def test_starts_at_zero(self):
        """Grid starts at z=0."""
        z = spatialGrid(0.05, 10, 5, 20)
        assert z[0] == 0.0

    def test_monotonically_increasing(self):
        """Grid depths are strictly increasing."""
        z = spatialGrid(0.05, 10, 5, 20)
        assert np.all(np.diff(z) > 0)

    def test_first_layer_thickness(self):
        """First layer thickness is zs/m."""
        zs = 0.05
        m = 10
        n = 5
        z = spatialGrid(zs, m, n, 20)
        dz0 = z[1] - z[0]
        np.testing.assert_allclose(dz0, zs / m, rtol=1e-10)

    def test_reaches_target_depth(self):
        """Grid extends to at least b skin depths."""
        zs = 0.05
        b = 20
        z = spatialGrid(zs, 10, 5, b)
        assert z[-1] >= zs * b

    def test_growth_ratio(self):
        """Layer thickness grows by factor (1+1/n)."""
        zs = 0.05
        n = 5
        z = spatialGrid(zs, 10, n, 20)
        dz = np.diff(z)
        ratios = dz[1:] / dz[:-1]
        expected = 1.0 + 1.0 / n
        np.testing.assert_allclose(ratios, expected, rtol=1e-10)


class TestInsertCustomLayerBoundaries:

    def test_no_op_without_custom_layers(self):
        z = spatialGrid(0.05, 10, 5, 20)
        result = insert_custom_layer_boundaries(z, [])
        assert result is z  # not just equal -- untouched, zero overhead

    def test_no_op_when_boundaries_already_present(self):
        """z_top=0 always coincides with node 0; a z_bottom that already
        matches an existing node needs no new node."""
        z = spatialGrid(0.05, 10, 5, 20)
        layer = DepthLayer(z_top=0.0, z_bottom=float(z[3]), rho=1100.0,
                           kc=0.001, chi=2.7)
        result = insert_custom_layer_boundaries(z, [layer])
        np.testing.assert_array_equal(result, z)

    def test_inserts_node_at_unmatched_boundary(self):
        """A boundary strictly between two nodes gets its own node,
        splitting that cell into two -- everything else is unchanged."""
        z = spatialGrid(0.05, 10, 5, 20)
        z_bottom = 0.5 * (z[0] + z[1])  # strictly inside the first cell
        layer = DepthLayer(z_top=0.0, z_bottom=z_bottom, rho=1100.0,
                           kc=0.001, chi=2.7)
        result = insert_custom_layer_boundaries(z, [layer])
        assert len(result) == len(z) + 1
        assert z_bottom in result
        assert np.all(np.diff(result) > 0)
        # Original nodes beyond the split are preserved exactly.
        np.testing.assert_array_equal(result[2:], z[1:])

    def test_inserts_both_boundaries_of_a_buried_layer(self):
        z = spatialGrid(0.05, 10, 5, 20)
        z_top = 0.5 * (z[2] + z[3])
        z_bottom = 0.5 * (z[3] + z[4])
        layer = DepthLayer(z_top=z_top, z_bottom=z_bottom, rho=1100.0,
                           kc=0.001, chi=2.7)
        result = insert_custom_layer_boundaries(z, [layer])
        assert len(result) == len(z) + 2
        assert z_top in result and z_bottom in result
        assert np.all(np.diff(result) > 0)

    def test_shared_boundary_between_adjacent_layers_inserted_once(self):
        z = spatialGrid(0.05, 10, 5, 20)
        z_mid = 0.5 * (z[0] + z[1])
        top_layer = DepthLayer(z_top=0.0, z_bottom=z_mid, rho=1100.0,
                               kc=0.001, chi=2.7)
        bottom_layer = DepthLayer(z_top=z_mid, z_bottom=float(z[1]),
                                  rho=1800.0, kc=0.003, chi=2.7)
        result = insert_custom_layer_boundaries(z, [top_layer, bottom_layer])
        assert len(result) == len(z) + 1  # one new node, not two
        assert np.sum(result == z_mid) == 1

    def test_near_duplicate_boundaries_merge_instead_of_creating_tiny_cell(self):
        """Two boundaries closer together than the merge tolerance must
        not create a vanishingly thin cell (which could force an
        arbitrarily small stable time step)."""
        z = spatialGrid(0.05, 10, 5, 20)
        z_mid = 0.5 * (z[0] + z[1])
        eps = 1e-9 * z[-1]  # far below the default 1e-6 * z[-1] tolerance
        layer_a = DepthLayer(z_top=0.0, z_bottom=z_mid, rho=1100.0,
                             kc=0.001, chi=2.7)
        layer_b = DepthLayer(z_top=z_mid + eps, z_bottom=float(z[1]),
                             rho=1800.0, kc=0.003, chi=2.7)
        result = insert_custom_layer_boundaries(z, [layer_a, layer_b])
        assert np.all(np.diff(result) > 1e-7 * z[-1])
        assert len(result) == len(z) + 1  # the two near-duplicates merged

    def test_boundary_beyond_grid_bottom_ignored(self):
        z = spatialGrid(0.05, 10, 5, 20)
        layer = DepthLayer(z_top=0.0, z_bottom=z[-1] + 1.0, rho=1100.0,
                           kc=0.001, chi=2.7)
        result = insert_custom_layer_boundaries(z, [layer])
        np.testing.assert_array_equal(result, z)

    def test_no_pathologically_tiny_cells_across_grid_resolutions(self):
        """A boundary landing a hair past a pre-existing background node
        (pure coincidence of the geometric grid spacing, not an
        intentional near-duplicate) must not leave that hair-thin gap
        unmerged: computeCFL's stable time step is limited by the
        *smallest* cell in the whole grid, so a single sliver cell can
        force a run to take orders of magnitude longer regardless of
        where the sliver sits. Sweeping the resolution parameter finds
        boundary/grid coincidences of this kind organically."""
        layer = DepthLayer(z_top=0.0, z_bottom=0.002, rho=1100.0,
                           kc=0.0004, chi=2.7)
        for m in (10, 20, 30, 50, 70, 100, 200, 400, 800):
            z = spatialGrid(0.045, m, 4, 20)
            result = insert_custom_layer_boundaries(z, [layer])
            dz = np.diff(result)
            dz_orig_min = np.diff(z).min()
            # No cell drops below ~10% of the ORIGINAL grid's own
            # smallest cell -- a modest, bounded CFL penalty regardless
            # of where the boundary happens to coincide with a node.
            assert dz.min() > 0.05 * dz_orig_min, (
                f"m={m}: sliver cell {dz.min():.3e} m vs original "
                f"minimum {dz_orig_min:.3e} m"
            )

    def test_result_always_sorted_and_unique(self):
        z = spatialGrid(0.05, 10, 5, 20)
        layers = [
            DepthLayer(z_top=0.0, z_bottom=0.5 * (z[0] + z[1]),
                      rho=1100.0, kc=0.001, chi=2.7),
            DepthLayer(z_top=0.5 * (z[4] + z[5]),
                      z_bottom=0.5 * (z[6] + z[7]),
                      rho=1800.0, kc=0.003, chi=2.7),
        ]
        result = insert_custom_layer_boundaries(z, layers)
        assert np.all(np.diff(result) > 0)
        assert len(np.unique(result)) == len(result)
