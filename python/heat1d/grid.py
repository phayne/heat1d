"""Spatial grid construction for heat1d.

Implements the non-uniform finite difference grid described in
Hayne et al. (2017), Appendix A2.2 (Eqs. A31-A33).
"""
import numpy as np


def skinDepth(P, kappa):
    """Calculate thermal skin depth.

    Parameters
    ----------
    P : float
        Period (e.g., diurnal, seasonal) [s]
    kappa : float
        Thermal diffusivity = k/(rho*cp) [m2.s-1]

    Returns
    -------
    float
        Thermal skin depth [m]
    """
    return np.sqrt(kappa * P / np.pi)


def spatialGrid(zs, m, n, b):
    """Calculate the spatial grid.

    The spatial grid is non-uniform, with layer thickness increasing
    geometrically downward: dz[i] = dz[0] * r^i where r = 1 + 1/n.

    Parameters
    ----------
    zs : float
        Thermal skin depth [m]
    m : int
        Number of layers in upper skin depth [default: 10, set in Config]
    n : int
        Layer increase with depth: dz[i] = dz[i-1]*(1+1/n) [default: 5, set in Config]
    b : int
        Number of skin depths to bottom layer [default: 20, set in Config]

    Returns
    -------
    np.ndarray
        Spatial node coordinates in meters.
    """
    dz0 = zs / m
    r = 1.0 + 1.0 / n
    zmax = zs * b
    # Compute number of layers analytically from the geometric series
    N = int(np.ceil(np.log(1 + zmax * (r - 1) / dz0) / np.log(r)))
    # Layer thicknesses: geometric progression
    dz = dz0 * r ** np.arange(N)
    # Node positions: cumulative sum with z[0] = 0
    z = np.zeros(N + 1)
    z[1:] = np.cumsum(dz)
    return z


def insert_custom_layer_boundaries(z, custom_layers, snap_fraction=0.1):
    """Insert grid nodes at every custom-layer boundary.

    :func:`spatialGrid` builds a grid from the *background* material's
    skin depth alone, with no knowledge of ``custom_layers``
    (see :mod:`heat1d.layers`). Left as-is, a layer whose ``z_top`` or
    ``z_bottom`` falls strictly between two grid nodes has its
    boundary *smeared* across whichever cell straddles it: that cell's
    conductivity is set from a single endpoint's value (see Eq A15 in
    :mod:`heat1d.solvers`), so the model ends up representing the
    layer's true thickness with a cell that is too thick or too thin,
    biasing its effective thermal resistance -- or, if the layer falls
    entirely inside one cell touching neither of its nodes,
    :func:`heat1d.layers.apply_custom_layers` misses it completely (no
    node satisfies its containment test).

    This function eliminates both failure modes by construction: it
    guarantees a grid node exists at every layer's ``z_top`` and
    ``z_bottom`` (splitting whichever cell used to straddle that depth
    into two), so every cell ends up either entirely inside or entirely
    outside each layer. Since ``DepthLayer`` properties are uniform
    within a layer, this makes the model's per-cell resistance exactly
    match the true material's, for a layer of any thickness or
    position -- see :func:`heat1d.diagnostics.layer_resolution_check`.

    A boundary landing at (or very close to) an existing node reuses it
    rather than inserting a near-duplicate. This matters beyond just
    avoiding redundant nodes: an explicit or Crank-Nicolson solver's
    stable time step is CFL-limited by the *smallest* cell in the whole
    grid (``computeCFL``, roughly ``dt <~ dz^2``), so even one
    accidental sliver cell -- e.g. a pre-existing background node
    landing a fraction of a percent past a layer boundary purely by
    coincidence of the geometric grid spacing -- can force a stable
    time step orders of magnitude smaller for the *entire* run, not
    just near that cell. "Close" is judged *locally*: a candidate node
    closer to its neighbor than ``snap_fraction`` times that
    neighboring cell's own width is merged into it, so the threshold
    scales naturally with the grid's geometric growth rather than using
    one absolute distance for both a sub-millimeter surface cell and a
    centimeter-scale cell near the bottom. When a requested boundary
    and an existing node are merged, the exact boundary depth wins
    (keeping the layer's edge precise) unless the existing node is
    itself another layer's boundary, in which case the earlier one is
    kept and the later one shifts by at most the local snap distance.
    Merging never moves a kept node by more than ``snap_fraction`` of
    its neighboring cell width, so it cannot itself create a large
    resistance error (see :func:`heat1d.diagnostics.layer_resolution_check`,
    which will still flag a layer whose thickness is itself smaller
    than this tolerance and so collapses to nothing).

    A boundary at exactly ``z=0`` (an ordinary surface-anchored layer)
    or at or beyond the grid's deepest node needs no insertion, since
    node 0 already sits there or no cell reaches that depth anyway.

    Parameters
    ----------
    z : np.ndarray
        Grid node depths [m], as returned by :func:`spatialGrid`.
    custom_layers : list of heat1d.layers.DepthLayer
    snap_fraction : float
        Merge a candidate node into its neighbor when closer than this
        fraction of the local (already-kept) cell width, rather than
        keep it as a separate, possibly very thin, cell.

    Returns
    -------
    np.ndarray
        Grid node depths [m], with layer boundaries included, sorted
        and deduplicated. Identical to ``z`` if ``custom_layers`` is
        empty or none of its boundaries fall strictly inside the grid.
    """
    if not custom_layers:
        return z

    z_max = z[-1]
    boundary_depths = sorted({
        float(zb)
        for layer in custom_layers
        for zb in (layer.z_top, layer.z_bottom)
        if 0.0 < zb < z_max
    })
    if not boundary_depths:
        return z

    combined = np.union1d(z, np.array(boundary_depths))
    is_boundary = np.isin(combined, np.array(boundary_depths))

    merged = [combined[0]]
    merged_is_boundary = [bool(is_boundary[0])]
    for i in range(1, len(combined)):
        # Local reference scale: the most recently kept cell width, or
        # (before any cell exists yet) the next gap in the combined
        # array, so a candidate immediately after z[0] is judged
        # against a sensible neighboring scale rather than nothing.
        if len(merged) >= 2:
            local_dz = merged[-1] - merged[-2]
        else:
            local_dz = combined[min(i + 1, len(combined) - 1)] - merged[-1]

        if combined[i] - merged[-1] < snap_fraction * local_dz:
            # Too close, relative to the local grid density, to resolve
            # as a separate cell; keep whichever is an exact layer
            # boundary so the layer's edge stays precise.
            if is_boundary[i] and not merged_is_boundary[-1]:
                merged[-1] = combined[i]
                merged_is_boundary[-1] = True
            # else: drop combined[i] (merge into the kept node)
        else:
            merged.append(combined[i])
            merged_is_boundary.append(bool(is_boundary[i]))

    return np.array(merged)
