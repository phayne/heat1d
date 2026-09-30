"""Energy-conservation diagnostics for heat1d.

Instruments the finite-difference scheme of Hayne et al. (2017),
Appendix A2 (Eqs. A13-A18, A23-A25, A30), to expose *where* — at which
depth and interface — the discrete model gains or loses energy, and how
that error accumulates over the full column.

Motivation
----------
The interior stencil (Eqs A15-A18) evaluates the conductive flux across
each inter-nodal segment using the conductivity at a *single* endpoint:

    F_j = k_j * (T[j+1] - T[j]) / dz_j        (Eq A15)

i.e. the shallower node's conductivity represents the whole segment.
This is a good approximation when k(z) varies slowly across a cell (the
smooth exponential profile of Eq A2 that this scheme was designed for),
but it is a poor one wherever a grid cell straddles a sharp property
contrast — most commonly a thin, low-conductivity custom layer (see
:mod:`heat1d.layers`) that the background-only grid (see
:func:`heat1d.grid.spatialGrid`, sized from ``planet.ks``/``rhos``/
``cp0`` alone, blind to any ``custom_layers``) fails to resolve.

Two more approximations compound this at the very top of the column:

* The surface boundary condition (Eqs A23-A25) estimates the flux at
  z=0 with a *three-point* one-sided derivative spanning nodes 0-2,
  assuming k is ~uniform over that span. The interior stencil instead
  feeds node 1 with the cruder *two-point* estimate ``F[0]`` above.
  These two flux estimates are only guaranteed to agree in the
  continuum limit; :func:`bc_flux_mismatch` quantifies their gap.
* Node 0 itself carries no explicit heat capacity in this scheme --
  ``surfTemp`` solves a quasi-steady balance, not a prognostic
  equation -- so it has no discrete storage term to check for
  conservation against; energy accounting below deliberately excludes
  it (see :func:`column_energy_budget`).

None of this is a bug in the ordinary (smoothly-varying, homogeneous or
gently layered) case the model was validated against. It becomes a
first-order effect for a thin, sharply-contrasting custom layer, which
is exactly the failure mode this module is built to localize.

What the flux-divergence probe actually finds
----------------------------------------------
Applied to a thin, low-conductivity custom layer thinner than the
grid's near-surface cell size, :func:`column_energy_budget` shows *no*
interior-node energy-conservation violation: ``residual`` sits at
floating-point noise at every node beyond the boundary-adjacent node 1,
whether or not a custom layer is present. The interior scheme is
exactly, locally conservative for whatever piecewise-constant material
properties it ends up using -- it is not "leaking" energy.

The reported inaccuracy instead comes from *which* properties it uses:
:func:`layer_resolution_check` shows that a grid cell straddling the
layer's boundary gets the wrong single-valued conductivity for its
entire span (the model has no way to represent a boundary that falls
between two nodes), which can misstate that cell's thermal resistance
by a large factor -- a *biased effective medium*, correctly integrated,
rather than a bookkeeping error. Use ``layer_resolution_check`` first
(no simulation needed) to check whether a layer is resolved at all, and
``column_energy_budget`` to confirm the interior scheme's own
conservation is intact regardless.

References
----------
Hayne, P. O. et al. (2017), JGR Planets 122, Appendix A.
"""
import numpy as np

from .layers import compute_default_properties
from .properties import heatCapacity, thermCond


def node_heat_capacity(planet, T, cp_model="polynomial"):
    """Heat capacity at every node for an arbitrary temperature array.

    Thin wrapper around :func:`heat1d.properties.heatCapacity`, named to
    pair with :func:`node_conductivity` -- both exist so a diagnostic can
    reconstruct T-dependent properties at a *saved* temperature state
    without relying on ``profile.k``/``profile.cp``, which only ever
    hold the values for the profile's *current* (typically final-step)
    temperature.

    Parameters
    ----------
    planet : object
        Planet object (as stored on ``model.planet``).
    T : np.ndarray
        Temperature [K].
    cp_model : str
        Heat capacity model, as in ``Configurator.cp_model``.

    Returns
    -------
    np.ndarray
        Heat capacity [J/kg/K], same shape as ``T``.
    """
    return heatCapacity(planet, np.asarray(T, dtype=float), model=cp_model)


def node_conductivity(profile, T):
    """Thermal conductivity at every node for an arbitrary temperature array.

    Reproduces :meth:`heat1d.profile.Profile.update_k` without mutating
    ``profile`` or requiring ``T`` to be the profile's *current*
    temperature -- so it can be evaluated against a saved ``T(t, z)``
    history from a completed run.

    Parameters
    ----------
    profile : heat1d.profile.Profile
        Supplies the (possibly per-node, if custom layers were used)
        contact conductivity ``kc`` and radiative parameter.
    T : np.ndarray
        Temperature [K], shape ``(N,)`` matching ``profile.z``.

    Returns
    -------
    np.ndarray
        Conductivity [W/m/K], shape ``(N,)``.
    """
    T = np.asarray(T, dtype=float)
    if profile._R350_array is not None:
        return profile.kc * (1.0 + profile._R350_array * T ** 3)
    return thermCond(profile.kc, T, profile.R350)


def conductive_flux(profile, T=None):
    """Conductive flux across every inter-nodal segment (Eq A15).

    ``F[j]`` is the flux between node ``j`` and node ``j+1``, evaluated
    with the shallower node's conductivity -- exactly the quantity the
    interior stencil (:func:`heat1d.solvers.solve_explicit` et al.) uses
    to advance the temperature field. ``F[-1]`` is also the flux the
    bottom boundary condition (Eq A30) enforces to equal ``planet.Qb``.

    Parameters
    ----------
    profile : heat1d.profile.Profile
    T : np.ndarray, optional
        Temperature [K], shape ``(N,)``. Defaults to ``profile.T``.

    Returns
    -------
    np.ndarray
        Flux [W/m^2], shape ``(N-1,)``, positive downward (increasing z).
    """
    T = profile.T if T is None else np.asarray(T, dtype=float)
    k = node_conductivity(profile, T)
    return k[:-1] * np.diff(T) / profile.dz


def surface_flux_3point(profile, T=None):
    """Conductive flux at z=0 from the boundary condition's 3-point stencil.

    Reconstructs the same quantity ``surfTemp`` (Eqs A24-A25) balances
    against absorbed and radiated flux, using nodes 0, 1, and 2. This is
    *not*, in general, equal to ``conductive_flux(profile, T)[0]`` (the
    2-point estimate the interior stencil actually feeds to node 1) --
    see :func:`bc_flux_mismatch`.

    Parameters
    ----------
    profile : heat1d.profile.Profile
    T : np.ndarray, optional
        Temperature [K], shape ``(N,)``. Defaults to ``profile.T``.

    Returns
    -------
    float
        Flux [W/m^2] at z=0.
    """
    T = profile.T if T is None else np.asarray(T, dtype=float)
    k0 = node_conductivity(profile, T[:1])[0]
    dTdz0 = (-3.0 * T[0] + 4.0 * T[1] - T[2]) / (2.0 * profile.dz[0])
    return k0 * dTdz0


def bc_flux_mismatch(profile, T=None):
    """Gap between the boundary's 3-point flux and the interior's 2-point flux.

    ``surface_flux_3point() - conductive_flux()[0]``, in W/m^2. Zero in
    the continuum limit for any smooth k(z); large and grid-dependent
    whenever k(z) changes sharply within the span the 3-point stencil
    covers (nodes 0-2) -- the diagnostic signature of an unresolved
    thin custom layer at the surface.

    Parameters
    ----------
    profile : heat1d.profile.Profile
    T : np.ndarray, optional
        Temperature [K], shape ``(N,)``. Defaults to ``profile.T``.

    Returns
    -------
    float
        Flux mismatch [W/m^2].
    """
    T = profile.T if T is None else np.asarray(T, dtype=float)
    F = conductive_flux(profile, T)
    return surface_flux_3point(profile, T) - F[0]


def flux_divergence(profile, T=None):
    """Conductive flux divergence at every interior node (Eq A16).

    Equal, node for node, to the bracketed term
    ``alpha*T[i-1] - (alpha+beta)*T[i] + beta*T[i+1]`` the interior
    solvers compute internally -- i.e. this is exactly the "flux
    divergence" the PDE (Eq A13, ``rho*cp*dT/dt = dF/dz``) attributes to
    conduction alone at each node, with no approximation beyond what
    the solver itself already makes.

    Parameters
    ----------
    profile : heat1d.profile.Profile
    T : np.ndarray, optional
        Temperature [K], shape ``(N,)``. Defaults to ``profile.T``.

    Returns
    -------
    np.ndarray
        ``dF/dz`` [W/m^3], shape ``(N-2,)``, for nodes ``1..N-2``.
    """
    T = profile.T if T is None else np.asarray(T, dtype=float)
    F = conductive_flux(profile, T)
    dz_cv = 0.5 * (profile.dz[1:] + profile.dz[:-1])
    return (F[1:] - F[:-1]) / dz_cv


def layer_probe(profile, T=None):
    """Snapshot of every quantity above, for a single temperature profile.

    Convenience wrapper bundling :func:`conductive_flux`,
    :func:`flux_divergence`, :func:`surface_flux_3point`, and
    :func:`bc_flux_mismatch` for one instant, e.g. for inspecting a
    live ``model.profile`` mid-run or a single saved output frame.

    Parameters
    ----------
    profile : heat1d.profile.Profile
    T : np.ndarray, optional
        Temperature [K], shape ``(N,)``. Defaults to ``profile.T``.

    Returns
    -------
    dict
        ``z`` (node depths, m), ``F`` (interface flux, W/m^2, length
        N-1), ``dFdz`` (interior flux divergence, W/m^3, length N-2),
        ``F_surf_3pt`` (float, W/m^2), ``bc_mismatch`` (float, W/m^2),
        ``Qb`` (float, W/m^2, ``F[-1]`` -- should equal
        ``profile.planet.Qb`` after ``botTemp`` has run).
    """
    T = profile.T if T is None else np.asarray(T, dtype=float)
    F = conductive_flux(profile, T)
    return {
        "z": profile.z,
        "F": F,
        "dFdz": flux_divergence(profile, T),
        "F_surf_3pt": surface_flux_3point(profile, T),
        "bc_mismatch": surface_flux_3point(profile, T) - F[0],
        "Qb": F[-1],
    }


def column_energy_budget(model):
    """Time-resolved conductive flux, flux divergence, and energy budget.

    Reconstructs, from a *completed* run's saved ``model.T(t, z)`` (no
    re-simulation, no modification to ``Model`` needed), everything
    :func:`layer_probe` reports at every output step, plus:

    * ``dTdt`` -- forward time derivative of the saved interior
      temperatures [K/s] (``(T[i+1]-T[i])/dt``, matching the update
      the explicit solver itself performs from state ``T[i]``), from
      which
    * ``residual`` -- ``rho*cp*dTdt - dFdz`` [W/m^3] at every interior
      node and output step. This is the discrete PDE residual: zero up
      to the *solver's own* temporal truncation error for a
      well-resolved column, and sharply nonzero at any node whose
      *spatial* discretization (grid cell straddling a property
      discontinuity) cannot represent the true flux divergence there.
    * ``column_storage_rate`` -- ``d/dt`` of the column's stored
      thermal energy, integrated over nodes ``1..N-2`` only (node 0 has
      no discrete heat capacity in this scheme; see module docstring).
    * ``column_residual`` -- ``(F[0] - F[-1]) - column_storage_rate``
      [W/m^2]: the column-integrated form of ``residual`` above (its
      node-centered control-volume sum), a stricter, time-resolved
      companion to :func:`heat1d.validation.check_energy_conservation`
      (which only compares the first and last output steps).
    * ``bc_mismatch`` -- time series of :func:`bc_flux_mismatch`: the
      part of any imbalance that originates at the surface-adjacent
      cell specifically, as opposed to the interior.

    Time derivatives use the actual sample spacing implied by
    ``model.lt`` (converted to seconds via ``model.planet.day``), not
    the solver's internal step size. ``dFdz`` at step ``i`` is
    evaluated from ``T[i]`` (the pre-step state), so for the explicit
    solver with every internal step saved (``config.solver="explicit"``,
    ``config.output_interval=None``) ``residual`` reproduces exactly
    what :func:`heat1d.solvers.solve_explicit` computed -- zero to
    floating-point precision for a well-resolved column, which is the
    cleanest configuration for isolating a *spatial* discretization
    defect (an unresolved custom layer) from time-integration error.
    Coarser output sampling, or the Crank-Nicolson/implicit solvers
    (whose sub-steps are not individually recorded), introduce their
    own temporal truncation error into ``residual`` alongside whatever
    spatial effect is being probed; that baseline is small and roughly
    uniform with depth, so a large, depth-*localized* residual still
    reliably points at the interface responsible. The final output
    step is dropped (no forward pair available there).

    Limitations
    -----------
    Node 1 (the shallowest interior node) always carries an intrinsic
    baseline residual, even for a smooth, homogeneous, well-resolved
    column: ``surfTemp`` solves the surface temperature for the *new*
    time level before the interior solver runs, so node 1's explicit
    update mixes the new ``T[0]`` with old ``T[1]``, ``T[2]`` -- a
    deliberate, documented trade-off of the operator-splitting scheme
    (see :mod:`heat1d.boundary`), not a defect. To isolate a genuine
    spatial-discretization problem (e.g. an unresolved custom layer)
    from this baseline, compare against a homogeneous run at the same
    grid resolution: the baseline is confined to node 1 and decays to
    numerical noise by node 2, whereas a real problem elevates the
    residual over multiple nodes and does not vanish with depth nearly
    as fast.

    ``Profile.update_properties`` skips recomputing ``k``/``cp`` when
    the temperature has changed by less than ``_prop_cache_threshold``
    (1 K by default) since the last update, so the *actual* solver used
    slightly stale properties that this function's fresh
    ``node_conductivity``/``node_heat_capacity`` evaluations do not
    reproduce exactly. For smooth media this adds a few W/m^3 of noise
    beyond node 1 (order 0.1-1% of the diurnal flux scale); set
    ``profile._prop_cache_threshold = 0.0`` before the run for an exact
    match when using this diagnostic.

    Assumes a flat surface with no PSR crater and no ground-heating
    term (``model.ground_heating``); ``Qb`` is taken as the constant
    ``model.planet.Qb``. Extending to those cases would need the
    corresponding extra source terms folded into ``column_residual``.

    Parameters
    ----------
    model : heat1d.model.Model
        A completed run (``model.run()`` already called).

    Returns
    -------
    dict
        ``t`` (s, length n-1), ``lt`` (hr, length n-1), ``z``,
        ``F`` (n-1, N-1), ``dFdz`` (n-1, N-2), ``dTdt`` (n-1, N-2),
        ``residual`` (n-1, N-2), ``bc_mismatch`` (n-1,),
        ``column_storage_rate`` (n-1,), ``column_residual`` (n-1,),
        ``Qb`` (n-1,, reconstructed top-of-bottom-cell flux), and
        ``Qb_expected`` (``model.planet.Qb``, float). Each series at
        index ``i`` describes the interval starting at output step
        ``i`` (i.e. is evaluated from ``model.T[i]``).
    """
    p = model.profile
    T_hist = model.T
    n_steps = T_hist.shape[0]
    if n_steps < 2:
        raise ValueError(
            "column_energy_budget needs at least 2 saved output steps "
            "to form a time derivative"
        )

    t_s = model.lt * (model.planet.day / 24.0)
    dt_fwd = t_s[1:] - t_s[:-1]  # length n-1, seconds spanning each interval

    dz_cv = 0.5 * (p.dz[1:] + p.dz[:-1])  # length N-2, interior control volumes

    n_out = n_steps - 1
    N = p.nlayers
    F = np.empty((n_out, N - 1))
    dFdz = np.empty((n_out, N - 2))
    bc_mismatch = np.empty(n_out)
    cp_hist = np.empty((n_out, N))

    for i in range(n_out):
        T_i = T_hist[i]  # pre-step state; matches what the solver used
        F[i] = conductive_flux(p, T_i)
        dFdz[i] = flux_divergence(p, T_i)
        bc_mismatch[i] = surface_flux_3point(p, T_i) - F[i, 0]
        # heatCapacity is T-dependent; profile.cp only ever holds the
        # value for whatever T the profile last had (typically the
        # final step, once the run has finished) -- it must be
        # recomputed at each saved T_i, not read as a static array.
        cp_hist[i] = node_heat_capacity(model.planet, T_i, p.config.cp_model)

    # Forward dT/dt at interior nodes from the SAVED output history --
    # matches the explicit solver's own update exactly when every
    # internal step is saved (see docstring).
    dTdt = (T_hist[1:, 1:-1] - T_hist[:-1, 1:-1]) / dt_fwd[:, None]

    rho_i = p.rho[1:-1]
    cp_i = cp_hist[:, 1:-1]
    residual = rho_i * cp_i * dTdt - dFdz

    column_storage_rate = np.sum(rho_i * cp_i * dz_cv * dTdt, axis=1)
    column_residual = (F[:, 0] - F[:, -1]) - column_storage_rate

    return {
        "t": t_s[:-1],
        "lt": model.lt[:-1],
        "z": p.z,
        "F": F,
        "dFdz": dFdz,
        "dTdt": dTdt,
        "residual": residual,
        "bc_mismatch": bc_mismatch,
        "column_storage_rate": column_storage_rate,
        "column_residual": column_residual,
        "Qb": F[:, -1],
        "Qb_expected": model.planet.Qb,
    }


def _true_kc_at(z, planet, custom_layers):
    """True (continuum) contact conductivity at arbitrary depths.

    Same precedence rule as :func:`heat1d.layers.apply_custom_layers`:
    starts from the background exponential profile (Eq A2) and
    overrides with the *last* custom layer (in list order) containing
    each depth, without regard to the model's grid.

    Parameters
    ----------
    z : np.ndarray
        Depths [m], arbitrary (not necessarily grid nodes).
    planet : object
        Planet object with ``ks``, ``kd``, ``rhos``, ``rhod``, ``H``.
    custom_layers : list of heat1d.layers.DepthLayer

    Returns
    -------
    np.ndarray
        Contact conductivity [W/m/K], same shape as ``z``.
    """
    kc, _ = compute_default_properties(z, planet)
    kc = np.array(kc, dtype=float, copy=True)
    for layer in custom_layers:
        kc[layer.contains(z)] = layer.kc
    return kc


def layer_resolution_check(profile, n_sub=200):
    """How well the model's grid resolves each custom layer.

    Analysis of :func:`heat1d.grid.spatialGrid`: the grid spacing near
    the surface, ``dz0 = skin_depth(planet.day, ks/(rhos*cp0)) / m``, is
    sized from the *background* material alone (``planet.ks``,
    ``rhos``, ``cp0``) -- it has no awareness of ``custom_layers``.
    ``apply_custom_layers`` (see :mod:`heat1d.layers`) then assigns a
    layer's properties only to whichever *discrete nodes* happen to
    land inside it; it does not insert nodes at the layer's boundaries
    or otherwise account for a boundary that falls between two nodes.

    Since the interior stencil (Eq A15) uses a single nodal
    conductivity value for the flux across an *entire* inter-nodal
    segment, any cell whose true material composition is split between
    a custom layer and its surroundings gets replaced by a
    homogeneous cell at the wrong (over- or under-) resistance. This
    function quantifies that mismatch directly: for every grid cell a
    layer overlaps, it compares the model's resistance
    (``dz / kc[shallow node]``, exactly what the solver uses) to the
    cell's *true* resistance (``integral of 1/kc_true(z) dz`` across
    the cell, from the actual layered material). No time-stepping is
    involved -- this is purely a property of the grid and the
    requested layers, and can be checked before running a model at
    all.

    A large ratio here (e.g. a 2mm, 25x-lower-conductivity surface
    veneer against the default grid's ~4.5mm first cell) explains a
    biased-but-locally-conservative diurnal temperature error even
    though :func:`column_energy_budget` shows no interior-node energy
    leak -- the scheme conserves energy exactly for the (wrong)
    effective medium it ends up representing.

    The per-cell *ratio* pinpoints which cell misrepresents the layer
    and by how much, but it is not, on its own, a reliable predictor of
    the net surface-temperature error: a thin cell that straddles the
    layer boundary can carry a huge ratio (its true resistance is small
    while the model rounds it entirely into the more resistive
    material) while contributing little in absolute terms, so the
    resulting temperature bias vs. grid resolution is *not* generally
    monotonic. ``total_resistance_error`` -- the sum of
    ``|R_model - R_true|`` [K*m^2/W] over every affected cell -- tracks
    the net bias more reliably and is the recommended single number for
    judging "is this layer resolved well enough".

    Parameters
    ----------
    profile : heat1d.profile.Profile
        Must have been constructed with ``custom_layers``.
    n_sub : int
        Sub-samples per grid cell used to numerically integrate the
        true resistance (trapezoidal rule).

    Returns
    -------
    list of dict
        One entry per custom layer, each with: ``label``, ``z_top``,
        ``z_bottom``, ``thickness`` [m], ``n_nodes_inside`` (grid nodes
        strictly within the layer), ``cells`` (list of
        ``{"cell": (j, z_j, z_j+1), "R_model", "R_true", "ratio"}``
        for every grid cell overlapping the layer, ``ratio`` =
        ``R_model / R_true`` -- 1.0 means exactly resolved, greater
        than 1 means the model over-insulates that cell, less than 1
        means it under-insulates it), ``max_ratio_deviation`` (the
        overlapping cell whose ratio deviates furthest from 1, as
        ``max(ratio, 1/ratio)``), ``total_resistance_error`` (sum of
        ``|R_model - R_true|`` [K*m^2/W] across all overlapping cells),
        and ``signed_resistance_error`` (signed sum: positive means the
        grid net over-insulates the layer, negative under-insulates
        it). A layer with no overlapping cell at all (completely
        invisible to the model -- possible for a thin *buried* layer
        no node's segment reaches) reports ``cells=[]`` and both error
        metrics as ``inf``.
    """
    z = profile.z
    dz = profile.dz
    kc = profile.kc
    planet = profile.planet
    custom_layers = profile.custom_layers or []

    results = []
    for layer in custom_layers:
        n_nodes_inside = int(np.sum(layer.contains(z)))

        cells = []
        for j in range(len(dz)):
            z0, z1 = z[j], z[j + 1]
            # Does this cell overlap the layer's depth range at all?
            if z1 <= layer.z_top or z0 >= layer.z_bottom:
                continue
            z_sub = np.linspace(z0, z1, n_sub + 1)
            kc_true = _true_kc_at(z_sub, planet, custom_layers)
            inv_kc = 1.0 / kc_true
            R_true = np.sum(0.5 * (inv_kc[1:] + inv_kc[:-1]) * np.diff(z_sub))
            R_model = dz[j] / kc[j]
            ratio = R_model / R_true
            cells.append({
                "cell": (j, z0, z1),
                "R_model": R_model,
                "R_true": R_true,
                "ratio": ratio,
            })

        if cells:
            worst = max(cells, key=lambda c: max(c["ratio"], 1.0 / c["ratio"]))
            max_ratio_deviation = max(worst["ratio"], 1.0 / worst["ratio"])
            total_resistance_error = sum(
                abs(c["R_model"] - c["R_true"]) for c in cells
            )
            signed_resistance_error = sum(
                c["R_model"] - c["R_true"] for c in cells
            )
        else:
            # No grid cell overlaps the layer at all: it is completely
            # invisible to the model (can happen for a thin *buried*
            # layer that no node's [z_j, z_j+1) segment reaches).
            max_ratio_deviation = np.inf
            total_resistance_error = np.inf
            signed_resistance_error = np.inf

        results.append({
            "label": layer.label,
            "z_top": layer.z_top,
            "z_bottom": layer.z_bottom,
            "thickness": layer.thickness,
            "n_nodes_inside": n_nodes_inside,
            "cells": cells,
            "max_ratio_deviation": max_ratio_deviation,
            "total_resistance_error": total_resistance_error,
            "signed_resistance_error": signed_resistance_error,
        })
    return results
