# Spatial Grid

The finite difference grid uses geometrically increasing layer thicknesses,
following Hayne et al. (2017), Appendix A2.2 (Eqs. A31--A33).

## Thermal Skin Depth

The *thermal skin depth* $z_s$ is the characteristic depth of penetration
of a periodic temperature wave:

$$
z_s = \sqrt{\frac{\kappa P}{\pi}}
$$

where $P$ is the forcing period (e.g., the diurnal period) and
$\kappa = K / (\rho c_p)$ is the thermal diffusivity.

For the Moon, z_s ≈ 4–7 cm depending on the H-parameter and
temperature-dependent thermal properties (Hayne et al., 2017). The grid is
constructed using surface-minimum properties (~3 cm), which ensures
adequate resolution near the surface where gradients are steepest.

## Grid Construction

The grid starts at $z = 0$ (the surface) with an initial layer thickness:

$$
\Delta z_0 = \frac{z_s}{m}
$$

where $m$ is the number of layers within the first skin depth (default: 10).

Layer thickness grows geometrically with depth:

$$
\Delta z_{i+1} = \Delta z_i \left(1 + \frac{1}{n}\right)
$$

where $n$ controls the growth rate (default: 5). Larger $n$ gives
more uniform layers; smaller $n$ gives faster growth.

The grid extends to a total depth of $b$ skin depths (default: 20), ensuring
that the bottom boundary is far enough below the surface that diurnal temperature
variations are negligible.

## Grid Parameters

| Parameter | Default | Description |
|---|---|---|
| $m$ | 10 | Layers per skin depth |
| $n$ | 5 | Growth factor (dz[i+1] = dz[i]*(1+1/n)) |
| $b$ | 20 | Total depth in skin depths |

For the Moon, these defaults produce approximately 45 layers extending to a
depth of about 1 m.

## Custom Layers and Grid Resolution

`custom_layers` (see `heat1d.layers.DepthLayer`) let you override
density and conductivity within a depth range — e.g. a thin, low-
thermal-inertia veneer on top of a coarser background. The grid above
is built entirely from the *background* material's properties
(`planet.ks`, `rhos`, `cp0`); it has no knowledge of any custom layer.
`apply_custom_layers` then assigns a layer's properties to whichever
*discrete grid nodes* happen to fall inside it — it does not insert
nodes at the layer's boundaries, and it does not blend properties for a
cell that straddles one.

Two consequences follow directly, both worth checking before trusting a
result:

- **A layer thinner than the local cell size is misrepresented, not
  rejected.** If only the surface node lands inside the layer, the
  model applies the layer's conductivity across the *entire* first
  cell — which can be several times the layer's true thickness,
  producing an artificially over- or under-insulated surface,
  depending on how the cell boundary happens to fall relative to the
  layer. This is not a violation of energy conservation (the scheme
  still conserves energy exactly for whatever effective medium it ends
  up using) — it is a biased *representation* of the medium, and the
  resulting surface-temperature error does not shrink monotonically
  with grid refinement, since it depends on where node positions
  happen to land relative to the layer boundary at each resolution.
- **A layer entirely between two nodes is invisible to the simulation
  regardless of how it changes `kc`/`rho`.** If neither node bounding a
  cell falls inside `[z_top, z_bottom)`, `apply_custom_layers` never
  touches that cell's properties at all.

`heat1d.diagnostics.layer_resolution_check(profile)` quantifies this
directly and requires no simulation: for every grid cell a custom layer
overlaps, it compares the model's assumed thermal resistance
(`dz / kc` at that cell, exactly what the solver uses) to the cell's
*true* resistance (integrating the actual layered material). A large
ratio, or a nonzero `total_resistance_error` relative to the layer's
own resistance, means the grid needs refining — raise `config.m` (and
possibly `config.n`) until the near-surface cell size is well below the
thinnest custom layer. `Profile` emits a `UserWarning` automatically
when this check finds a layer resolved to worse than 50% of its own
resistance.

`heat1d.diagnostics.column_energy_budget(model)` complements this with
a time-resolved, per-layer flux-divergence probe on a completed run,
useful for confirming that a poorly-resolved layer is producing a
representation bias rather than an actual solver defect (see that
module's docstring for the full picture, including a known, unrelated
baseline residual at the node nearest the surface).
