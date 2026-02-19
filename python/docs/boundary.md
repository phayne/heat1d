# Boundary Conditions

The thermal model requires boundary conditions at the surface (top) and at depth
(bottom). These follow Hayne et al. (2017), Appendix A1.1 and A2.1
(Eqs. A7--A12, A21--A30).

## Surface Boundary Condition

The surface energy balance determines the surface temperature $T_s$:

$$
\varepsilon \sigma T_s^4 = Q_s + K \left. \frac{\partial T}{\partial z} \right|_{z=0}
$$

where:

- $\varepsilon$ is the infrared emissivity (0.95 for the Moon)
- $\sigma$ is the Stefan-Boltzmann constant
- $Q_s$ is the absorbed solar flux
- The last term is the conductive heat flux from below

The absorbed solar flux depends on the solar incidence angle $\theta$:

$$
Q_s = \frac{S_0}{r^2} (1 - A(\theta)) \cos\theta
$$

where $S_0$ is the solar constant at 1 AU, $r$ is the heliocentric
distance in AU, and $A(\theta)$ is the incidence angle-dependent albedo.

### Angle-Dependent Albedo

The albedo model follows Keihm (1984) and Vasavada et al. (2012), as used in
Hayne et al. (2017), Eq. A8:

$$
A(\theta) = A_0 + a \left(\frac{\theta}{\pi/4}\right)^3 +
b \left(\frac{\theta}{\pi/2}\right)^8
$$

where $A_0$ is the normal bolometric Bond albedo (0.12 for highland,
0.06 for mare; Hayne et al. 2017), $a = 0.06$, and $b = 0.25$
for the Moon. The angle-dependent enhancement causes the effective hemispheric
Bond albedo to exceed $A_0$.

### Newton's Method

The surface energy balance is a nonlinear equation in $T_s$ and is solved
iteratively using Newton's root-finding method. The iteration:

$$
T_s^{(k+1)} = T_s^{(k)} - \frac{f(T_s^{(k)})}{f'(T_s^{(k)})}
$$

converges when |ΔT| < ε (default: 0.1 K). A maximum
iteration count (default: 100) prevents infinite loops.

The surface temperature gradient is approximated using a second-order forward
difference:

$$
\left. \frac{\partial T}{\partial z} \right|_{z=0} \approx
\frac{-3 T_0 + 4 T_1 - T_2}{2 \Delta z_0}
$$

## Bottom Boundary Condition

The bottom boundary is set by the interior heat flux $Q_b$:

$$
T_N = T_{N-1} + \frac{Q_b}{K_{N-1}} \Delta z_{N-1}
$$

For the Moon, $Q_b = 0.018$ W m⁻² (Langseth et al., 1976).

## PSR Crater Illumination

Permanently shadowed regions (PSRs) on the floors of polar craters receive no
direct sunlight, but they are heated by sunlight scattered off the illuminated
crater walls and by thermal infrared radiation re-emitted by those walls. `heat1d`
models this using the analytical framework of
[Ingersoll & Svitek (1992)](https://doi.org/10.1016/0019-1035(92)90017-2)
for spherical-cap (bowl-shaped) craters.

### Crater Geometry

A bowl-shaped crater is parameterized by its depth-to-diameter ratio $d/D$.
The key derived quantity is the *surface area ratio* $f$, which measures the
fraction of the hemisphere subtended by the crater opening as seen from the
floor:

$$
f = \frac{4 (d/D)^2}{1 + 4 (d/D)^2}
$$

The corresponding crater half-angle $\beta$ (the angle from the crater axis to
the rim, as seen from the center of curvature) is:

$$
\beta = \arccos(1 - 2f)
$$

Typical simple lunar craters have d/D ≈ 0.1–0.2, giving
f ≈ 0.04–0.14 and β ≈ 23°–44°.

### PSR Viability Condition

A permanently shadowed region can exist only when the Sun never rises above the
crater rim as seen from the floor.  The maximum solar elevation angle at
latitude $\phi$ on a body with obliquity $\epsilon$ is:

$$
e_{0,\max} = \frac{\pi}{2} - \lvert\phi\rvert + \epsilon
$$

A PSR exists when:

$$
e_{0,\max} < \beta
$$

This sets a **minimum depth-to-diameter ratio** for PSR formation at a given
latitude. Inverting the geometry relations gives:

$$
\left(\frac{d}{D}\right)_{\min} = \frac{1}{2} \sqrt{\frac{1 - \cos e_{0,\max}}{1 + \cos e_{0,\max}}}
$$

For the Moon (ε ≈ 1.54°):

| Latitude | $e_{0,\max}$ | Min $d/D$ | Min $\beta$ |
|----------|---------------|-----------|-------------|
| 90° (pole) | 1.5° | 0.007 | 1.5° |
| 85° | 6.5° | 0.029 | 6.5° |
| 80° | 11.5° | 0.051 | 11.5° |
| 70° | 21.5° | 0.097 | 21.5° |

At lower latitudes the required crater depth increases rapidly, and below
~70° latitude even the deepest simple craters (d/D ≈ 0.2) cannot
sustain a PSR on the Moon.

If PSR mode is enabled at a latitude where the viability condition is not met,
the model issues a warning and runs anyway — the crater floor will receive
direct sunlight at high solar elevations, but the cavity scattering model still
applies.

### Cavity Radiation Trapping

The concave crater geometry traps infrared radiation: thermal photons emitted
by the floor are partially intercepted and re-absorbed by the walls, and vice
versa.  This increases the effective emissivity beyond the flat-surface value
(Ingersoll & Svitek 1992, Eq. 7):

$$
\varepsilon_{\text{eff}} = \frac{\varepsilon}{1 - (1 - \varepsilon) f}
$$

For $\varepsilon = 0.95$ and $f = 0.14$ ($d/D = 0.2$), the effective
emissivity increases to ~0.957.  The effect is modest for typical lunar
emissivities but becomes significant for lower-emissivity surfaces.

In the surface energy balance, $\varepsilon_{\text{eff}}$ replaces
$\varepsilon$ on the left-hand side:

$$
\varepsilon_{\text{eff}} \sigma T_s^4 = Q_{\text{psr}} + K \left. \frac{\partial T}{\partial z} \right|_{z=0}
$$

### Absorbed Flux at the Crater Floor

The absorbed solar flux at the PSR floor comes entirely from scattering and
re-emission by the sunlit crater walls (Ingersoll & Svitek 1992, Eq. 8):

$$
Q_{\text{psr}} = \frac{S_0}{r^2} \sin e_0 \cdot \frac{f (1 - A)}{1 - A f} \left[\varepsilon + A (1 - f)\right]
$$

where $e_0$ is the solar elevation angle (clamped to ≥ 0), $A$ is the
Bond albedo, and $\varepsilon$ is the flat-surface emissivity.

The three factors have clear physical meanings:

- $S_0 \sin e_0 / r^2$: solar flux intercepted by the crater opening (projected area)
- $f (1 - A) / (1 - Af)$: fraction absorbed by the cavity after multiple wall reflections
- $\varepsilon + A(1 - f)$: partition between thermal re-emission ($\varepsilon$) reaching the floor and reflected sunlight escaping through the opening ($A(1-f)$)

At night ($\sin e_0 = 0$), $Q_{\text{psr}} = 0$ and the crater floor cools
radiatively through the opening, governed by the effective emissivity.

### Solar Elevation Angle

The PSR flux model requires the solar elevation angle $e_0$ rather than the
incidence angle $\theta$ used for flat surfaces. For a horizontal surface at
latitude $\phi$ with solar declination $\delta$ and hour angle $h$:

$$
\sin e_0 = \cos\theta = \sin\phi \sin\delta + \cos\phi \cos\delta \cos h
$$

This uses the same orbital mechanics as the flat-surface model (see
[Theory](theory.md)). The PSR model clamps sin e₀ ≥ 0 so that negative
values (Sun below the horizon) produce zero flux rather than unphysical
negative flux.

### Integration with SPICE Ephemerides

When an external flux time series is provided (e.g., from JPL Horizons/SPICE
via `--use-spice`), the precomputed direct flux replaces the analytical
illumination model during the output phase. In the current implementation, PSR
mode and external flux series are **mutually exclusive**: the PSR cavity
scattering formula requires the solar elevation angle, which is computed from
the analytical orbit model rather than extracted from a precomputed flux value.

For SPICE-driven PSR simulations, the recommended workflow is:

1. Use the analytical orbit model with PSR mode for equilibration and output
2. Verify that the analytical ephemeris adequately reproduces the SPICE geometry
   at the latitude of interest

Future versions may support extracting solar elevation from SPICE geometry
kernels directly.

### Solver Compatibility

PSR mode is compatible with all three time-stepping solvers (explicit, implicit,
Crank-Nicolson) and with adaptive timestepping. The **Fourier-matrix solver is
not supported** for PSR because it requires the surface flux to be expressible
as a function of the surface temperature alone, whereas the PSR flux depends
nonlinearly on the solar elevation angle through the cavity scattering formula.
When the Fourier-matrix solver is selected with PSR mode enabled, the model
automatically falls back to the implicit solver.

### Usage

**Python API:**

```python
from heat1d import Model
import planets
import numpy as np

model = Model(planet=planets.Moon, lat=np.deg2rad(85), ndays=1, psr_d_D=0.2)
model.run()
print(f"PSR floor Tmax: {model.T[:, 0].max():.1f} K")
```

**CLI:**

```bash
heat1d --lat 85 --psr-d-D 0.2 --ndays 1
```

### Reference

Ingersoll, A. P., Svitek, T., & Murray, B. C. (1992). Stability of polar
frosts in spherical bowl-shaped craters on the Moon, Mercury, and Mars.
*Icarus*, 100, 40--47.
[doi:10.1016/0019-1035(92)90017-2](https://doi.org/10.1016/0019-1035(92)90017-2)
