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
(`planet.ks`, `rhos`, `cp0`); it has no knowledge of any custom layer
on its own.

That used to matter, because `apply_custom_layers` assigns a layer's
properties to whichever *discrete grid nodes* happen to fall inside it
— it does not, by itself, insert nodes at the layer's boundaries or
blend properties for a cell that straddles one. A layer thinner than
the local cell size would be misrepresented (its conductivity applied
across a whole cell several times its true thickness, biasing the
surface's effective resistance non-monotonically with grid
refinement), and a layer landing entirely between two nodes would be
invisible to the simulation regardless of its property contrast.

`Profile` now closes this gap automatically: before deriving grid
spacing, cell coefficients, or default properties, it calls
`heat1d.grid.insert_custom_layer_boundaries` to guarantee a node at
every custom layer's `z_top` and `z_bottom` (reusing an existing node,
or merging two nearly-coincident boundaries, when they already fall
close enough together). Every cell then lies entirely inside or
entirely outside each layer, so `apply_custom_layers`'s per-node
assignment gives every cell the *exact* resistance of the material it
actually contains — for a layer of any thickness or position, at any
grid coarseness (`config.m`/`config.n`/`config.b`). This does mean a
custom layer adds one or two extra nodes to the grid; it is otherwise
a no-op when no custom layers are given.

The one case this cannot fix is a layer thinner than its own merge
tolerance (a small fraction of the total grid depth) — it gets snapped
away as a near-duplicate node rather than resolved. `Profile` emits a
`UserWarning` automatically should this (or any other resolution
shortfall) occur.

`heat1d.diagnostics.layer_resolution_check(profile)` is the check this
warning is built on, and can be run directly with no simulation needed:
for every grid cell a custom layer overlaps, it compares the model's
assumed thermal resistance (`dz / kc` at that cell, exactly what the
solver uses) to the cell's *true* resistance (integrating the actual
layered material) — reporting `ratio = 1.0` and
`total_resistance_error ≈ 0` in the ordinary case.
`heat1d.diagnostics.column_energy_budget(model)` complements this with
a time-resolved, per-layer flux-divergence probe on a completed run,
confirming the interior scheme's own energy conservation is intact
independent of grid resolution (see that module's docstring for the
full picture, including a known, unrelated baseline residual at the
node nearest the surface).

**What this fix does not remove**: an exactly-represented layer still
needs *some* surrounding resolution to capture the diurnal transient
across a sharp property contrast accurately — this is ordinary
finite-difference truncation error, present (at a much smaller scale)
even for a smooth homogeneous profile, not something specific to custom
layers. For a 2mm, 25x-lower-conductivity surface veneer on Europa, the
default `m=10` still gives surface temperatures several K off from a
converged answer, improving roughly monotonically and converging to
<0.1K by `m~100`. Where the pre-fix behavior needed `m~400-800` *and*
converged non-monotonically (grid refinement could occasionally make
things worse, depending on where a node happened to land relative to
the layer boundary), the fixed grid converges smoothly and an order of
magnitude faster — but a very thin, sharply-contrasting layer still
benefits from checking a couple of `m` values against each other before
trusting the result.
