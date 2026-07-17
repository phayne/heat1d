# Sloped Surfaces

`heat1d` can model temperatures on a tilted planar surface — an isolated
slope surrounded by flat terrain — given the slope angle and azimuth.
The model accounts for:

1. **Direct insolation** on the tilted surface (Braun & Mitchell 1983
   geometry),
2. **Self-shadowing** — no direct flux when the sun is behind the slope
   or below the flat horizon,
3. **Ground heating** (optional, on by default) — indirect flux from
   the surrounding flat terrain: thermal emission plus reflected
   sunlight, computed from an automatic flat-surface companion run.

## Conventions

- **Slope** $s$: dip from horizontal, $0° \le s \le 90°$.
- **Azimuth** $\gamma$: the direction the tilted surface *faces* (the
  downslope direction), measured clockwise from north:
  0° = N, 90° = E, 180° = S, 270° = W. This matches the JPL Horizons
  `AZ` convention.
- User-facing interfaces (CLI, YAML, GUI, C config) use **degrees**;
  the Python `Model` API uses **radians** (like `lat`).

## Direct insolation geometry

With solar zenith angle $z$ and solar azimuth $A_\odot$ (computed from
latitude, declination, and hour angle), the cosine of the local
incidence angle on the tilted surface is

$$
\mu \;=\; \cos z \cos s \;+\; \sin z \sin s \cos(A_\odot - \gamma),
$$

subject to two shadowing conditions:

- $\mu < 0$: the sun is behind the slope (self-shadowing) → $\mu = 0$,
- $\cos z \le 0$: the sun is below the flat horizon → $\mu = 0$
  (a steep slope can satisfy $\mu > 0$ with the sun below the horizon;
  the surrounding terrain blocks it).

The absorbed direct flux is

$$
Q_\mathrm{dir} = \bigl[1 - A(i_\mathrm{loc})\bigr]\,
\frac{S}{(r/r_\mathrm{AU})^{2}}\, \mu,
\qquad i_\mathrm{loc} = \arccos\mu,
$$

where the angle-dependent albedo $A(i)$ (Keihm 1984; Vasavada et
al. 2012) is evaluated at the **local** incidence angle on the tilted
surface, since it is a photometric property of the regolith relative to
its own surface normal.

Note that on an airless body a slope can see the sun appear or vanish
**discontinuously** (e.g. an east-facing slope at sunrise receives
$\mu = \sin s$ of the full solar beam the moment the sun crests the
horizon). This is physical, but it produces localized Gibbs ringing in
the spectral `fourier-matrix` output solver near the terminator; a
time-stepping output solver is recommended for east/west-facing slopes
(the model warns when this applies). Pole- and equator-facing slopes
have continuous flux and are unaffected.

## View factors

An infinite tilted plane sees the surrounding flat ground with view
factor

$$
F_\mathrm{ground} = \sin^2(s/2), \qquad
F_\mathrm{sky} = \cos^2(s/2),
$$

Sky emission is negligible for airless bodies, so the sky view factor
does not contribute a flux term.

## Ground heating (indirect flux)

When `ground_heating` is enabled (the default for sloped runs), the
model first runs a **flat companion model** at the same latitude,
longitude, and configuration, extracts one equilibrated diurnal cycle
of its surface temperature $T_\mathrm{flat}(t)$, and analytically
recomputes the flat-terrain illumination. The indirect flux added to
the slope is the periodic series

$$
Q_\mathrm{ind}(t) = \sin^2(s/2)\,\Bigl[
\varepsilon^2 \sigma T_\mathrm{flat}^4(t)
\;+\; (1 - A_0)\, A(\theta(t))\, \frac{S}{(r/r_\mathrm{AU})^2}\,
\max(\cos\theta(t),\,0)
\Bigr]
$$

- **Thermal term**: the flat-ground radiosity
  $\varepsilon \sigma T^4$ is absorbed by the slope with infrared
  absorptivity $\varepsilon$ (Kirchhoff's law), giving the factor
  $\varepsilon^2$.
- **Reflected-solar term**: flat ground reflects
  $A(\theta) F_\mathrm{inc}$ of the incident sunlight (at solar zenith
  angle $\theta$); the slope absorbs the fraction $1 - A_0$ of this
  hemispherically diffuse illumination using the normal-incidence
  albedo $A_0$ (the grazing-incidence boost in $A(i)$ is a
  directional-beam effect).

$Q_\mathrm{ind}(t)$ is included during **both** equilibration and
output, indexed periodically by time-of-day. The same formulation is
implemented in the Python and C backends.

### Approximations

- **Infinite flat plane**: the surroundings are horizontal and extend
  to the horizon; local topography beyond the single slope is not
  modeled.
- **One-way coupling**: the slope does not heat the flat terrain back
  (second order in $\sin^2(s/2)$).
- **Single-cycle periodicity**: one equilibrated diurnal cycle is
  reused for all days — exact for circular orbits, and the same
  approximation the fourier-matrix equilibration makes for the direct
  flux. For eccentric orbits the day-to-day drift of $r$ and
  declination is not captured in $Q_\mathrm{ind}$.
- **Isotropic scattering**: reflected sunlight is treated as
  Lambertian.
- **Eclipses** (Horizons mode) reduce the direct beam only; the
  indirect flux is not eclipsed (a second-order effect).

## Usage

Command line:

```bash
# 25-degree slope facing south at latitude 30 N
heat1d --lat 30 --slope 25 --slope-az 180

# Disable the indirect terrain flux
heat1d --lat 30 --slope 25 --slope-az 180 --no-ground-heating

# With the C backend
heat1d --backend c --lat 30 --slope 25 --slope-az 180

# With JPL Horizons illumination (uses the queried solar azimuth/elevation)
heat1d --use-spice --lat 0 --slope 20 --slope-az 90 \
       --start-time "2024-06-01 12:00"
```

YAML (top-level keys, next to `latitude`):

```yaml
slope: 25.0            # degrees
slope_azimuth: 180.0   # degrees clockwise from N
ground_heating: true
```

Python API (radians):

```python
import numpy as np
from heat1d import Model, planets

m = Model(planet=planets.Moon, lat=np.deg2rad(30),
          slope=np.deg2rad(25), slope_az=np.deg2rad(180))
m.run()
```

The building blocks are exposed in `heat1d.terrain`
(`slope_incidence_cos`, `ground_view_factor`, `direct_slope_flux`,
`indirect_flux_series`) and `heat1d.orbits.solarAzimuth`.

Notes:

- `--slope` is mutually exclusive with `--psr-d-D` (different surface
  geometries).
- `--flux-file` cannot be combined with `--slope`: an external flux
  file is assumed to already contain slope-projected direct flux.
  Generate one with `generate-flux --slope <deg> --slope-az <deg>`.

## References

- Braun, J. E., & Mitchell, J. C. (1983). Solar geometry for fixed and
  tracking surfaces. *Solar Energy*, 31(5), 439–444.
- Aharonson, O., & Schorghofer, N. (2006). Subsurface ice on Mars with
  rough topography. *JGR*, 111, E11007.
- Keihm, S. J. (1984). Interpretation of the lunar microwave
  brightness temperature spectrum. *Icarus*, 60, 568–589.
- Vasavada, A. R., et al. (2012). Lunar equatorial surface temperatures
  and regolith properties from the Diviner Lunar Radiometer
  Experiment. *JGR*, 117, E00H18.
- Hayne, P. O., et al. (2017). Global regolith thermophysical
  properties of the Moon from the Diviner Lunar Radiometer Experiment.
  *JGR Planets*, 122, 2371–2400.
