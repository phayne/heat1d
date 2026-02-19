# heat1d

**One-dimensional thermal model for planetary science applications**

[![PyPI version](https://img.shields.io/pypi/v/heat1d.svg)](https://pypi.python.org/pypi/heat1d)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Documentation](https://readthedocs.org/projects/heat1d/badge/?version=latest)](https://heat1d.readthedocs.io)

`heat1d` solves the 1-D heat equation in planetary regolith with depth- and temperature-dependent thermophysical properties, following the formulation of [Hayne et al. (2017)](https://doi.org/10.1002/2017JE005387). Implementations are provided in Python, C, and MATLAB.

![Lunar Equatorial Surface Temperature](docs/images/diurnal_surface_temperature.png)

## Features

- **Four numerical solvers**: Explicit (Forward Euler), Implicit (Backward Euler), Crank-Nicolson, and Fourier-matrix (frequency domain)
- **Temperature-dependent properties**: thermal conductivity with radiative $T^3$ term, polynomial heat capacity
- **Depth-dependent profiles**: exponential density/conductivity transition with custom layer support
- **Angle-dependent albedo** following Keihm (1984)
- **Non-uniform spatial grid** with geometrically increasing layer spacing
- **Graphical user interface** (PySide6) with interactive parameter panel, run manager, and 5 plot types
- **JPL Horizons/SPICE integration** for precise ephemeris-driven illumination
- **Eclipse modeling** for satellite bodies (e.g., Moon in Earth's shadow)
- **PSR crater modeling**: bowl-shaped permanently shadowed regions ([Ingersoll & Svitek, 1992](docs/Ingersoll-Svitek_bowl-shaped-craters-frost_Icarus_1992.pdf))
- **YAML configuration** system shared between Python and C implementations
- **Validation suite** against Apollo heat flow data and Diviner radiometer observations
- **C implementation** for the highest computational performance (same YAML config, all solvers, 24 validation tests)

## Quick Start

### Installation

```bash
pip install heat1d
```

For the graphical user interface:

```bash
pip install "heat1d[gui]"
```

### Command Line

```bash
# Default: Moon equator, 1 diurnal cycle
heat1d

# Moon at 45 N with implicit solver, 2 output days
heat1d --lat 45 --solver implicit --ndays 2

# Load a YAML config file, override latitude
heat1d moon.yaml --lat 30

# Use JPL Horizons for real ephemeris
heat1d --use-spice --lat 0 --lon 0 --start-time "2024-06-15 12:00" --ndays 1

# Run the validation suite
heat1d --validate

# Launch the GUI
heat1d-gui
```

### Python API

```python
import numpy as np
import planets
from heat1d import Model, Configurator

config = Configurator(solver="crank-nicolson")
model = Model(planet=planets.Moon, lat=np.deg2rad(30), ndays=1, config=config)
model.run()

# Surface temperature vs local time
print(f"Peak T: {model.T[:, 0].max():.1f} K")
print(f"Min  T: {model.T[:, 0].min():.1f} K")

# Access depth grid and full T(time, depth) array
print(f"Grid: {model.N_z} layers to {model.profile.z[-1]:.2f} m depth")
```

## Example Output

### Temperature vs. Depth and Local Time

The thermal wave penetrates only a few centimeters into lunar regolith. Below the thermal skin depth (~7 cm, accounting for temperature-dependent conductivity and heat capacity), diurnal variations are strongly damped.

![Depth Heatmap](docs/images/depth_heatmap.png)

### Multi-Latitude Comparison

Peak daytime temperatures decrease with latitude due to lower solar incidence angles, while nighttime temperatures converge as all surfaces radiate to the same cold sky.

![Multi-Latitude](docs/images/multi_latitude.png)

## Solvers

`heat1d` implements four numerical methods for the 1-D heat equation. Three are finite-difference time-stepping schemes; the fourth solves in the frequency domain:

| Solver | Method | Stability | Time Accuracy | Typical Steps/Day | Relative Speed |
|---|---|---|---|---|---|
| **Explicit** | Forward Euler | CFL-limited | $O(\Delta t)$ | ~830 | 1x |
| **Implicit** | Backward Euler + TDMA | Unconditional | $O(\Delta t)$ | ~24 | ~35x |
| **Crank-Nicolson** | Semi-implicit + TDMA | Unconditional | $O(\Delta t^2)$ | ~24 | ~35x |
| **Fourier-matrix** | Frequency domain | N/A (periodic) | Spectral | N/A | ~1000x |

![Solver Comparison](docs/images/solver_comparison.png)

### Time-Stepping Solvers

The **explicit scheme** (Hayne et al. 2017, Eq. A17) is a straightforward forward Euler discretization. It updates each grid node from its neighbors at the current time step, requiring no linear algebra. However, the Courant-Friedrichs-Lewy (CFL) stability condition ($\Delta t \le \frac{\Delta z^2}{2\kappa}$) limits the time step to ~3000 seconds for typical lunar parameters, resulting in ~830 steps per synodic day.

The **implicit** and **Crank-Nicolson** schemes evaluate spatial derivatives at the new time level (fully or partially), producing a tridiagonal linear system at each step. This system is solved in O(N) time using the Thomas algorithm (TDMA). Because they are unconditionally stable, both schemes can take time steps ~35x larger than the explicit scheme. The Crank-Nicolson scheme averages explicit and implicit contributions, achieving second-order accuracy in time at the same computational cost as the first-order implicit scheme.

Both the implicit and Crank-Nicolson solvers support **adaptive time-stepping** via Richardson extrapolation (step doubling), automatically adjusting the time step to maintain a user-specified temperature accuracy.

### Fourier-Matrix Solver

The **Fourier-matrix solver** eliminates time-stepping entirely by solving the periodic steady-state in the frequency domain. It decomposes the diurnal surface flux into Fourier harmonics and propagates each frequency through the subsurface using 2×2 transmission matrices. Nonlinear surface radiation is handled via Newton iteration in the time domain, using a circulant admittance matrix constructed from the frequency-domain impedance.

An outer iteration loop captures the **solid-state greenhouse effect** (thermal pumping): the nonlinear $T^3$ dependence of thermal conductivity produces a net downward heat flux that elevates subsurface temperatures. The solver computes the exact time-averaged rectification flux and adjusts the equilibrium temperature profile accordingly.

This approach is ~1000× faster than time-stepping because it solves the periodic steady state directly, without equilibration orbits. A complete lunar diurnal cycle is solved in ~100 ms. It is also the **default equilibration solver** -- even when using time-stepping methods for output, the Fourier solver initializes the temperature profile to the converged periodic state, eliminating multi-orbit spin-up.

The Fourier-matrix approach builds on the harmonic decomposition of periodic lunar thermal models introduced by [Linsky (1966)](https://doi.org/10.1016/0019-1035(66)90075-3), the thermal quadrupole (transfer matrix) formalism of [Pipes (1957)](https://doi.org/10.1016/0016-0032(57)90927-6) and [Maillet et al. (2000)](https://doi.org/10.1002/9780470694374), and the harmonic balance technique for nonlinear periodic circuits ([Kundert & Sangiovanni-Vincentelli, 1986](https://doi.org/10.1109/TCAD.1986.1270223)). The solid-state greenhouse (thermal pumping) correction follows the analysis of [Linsky (1966)](https://doi.org/10.1016/0019-1035(66)90075-3), Section V. For detailed equations and derivations, see the [Numerical Methods](python/docs/numerical.md) documentation.

## PSR Crater Modeling

`heat1d` can model temperatures on the floors of **permanently shadowed regions** (PSRs) in bowl-shaped polar craters, following the analytical framework of [Ingersoll & Svitek (1992)](https://doi.org/10.1016/0019-1035(92)90017-2). The crater floor receives no direct sunlight but is heated by scattered sunlight and thermal infrared re-emitted by the illuminated crater walls.

The crater geometry is parameterized by the depth-to-diameter ratio $d/D$, which determines the surface area ratio $f = 4(d/D)^2 / [1 + 4(d/D)^2]$ and the crater half-angle $\beta = \arccos(1 - 2f)$. A PSR exists when the maximum solar elevation angle at the given latitude is less than $\beta$, setting a minimum $d/D$ that increases with distance from the pole. For the Moon ($\epsilon \approx 1.54$°), even shallow craters ($d/D \approx 0.01$) can host PSRs at the poles, while at 70° latitude only the deepest simple craters ($d/D \gtrsim 0.1$) qualify.

The model accounts for cavity radiation trapping (enhanced effective emissivity) and computes the absorbed flux from wall scattering. PSR mode works with all time-stepping solvers and adaptive timestepping, but not with the Fourier-matrix solver (which automatically falls back to implicit).

```bash
# PSR crater floor at 85°N with d/D = 0.2
heat1d --lat 85 --psr-d-D 0.2 --ndays 1
```

```python
model = Model(planet=planets.Moon, lat=np.deg2rad(85), ndays=1, psr_d_D=0.2)
model.run()
```

For the full derivation and formulas, see the [Boundary Conditions](python/docs/boundary.md) documentation.

## Graphical User Interface

The GUI provides interactive control over all model parameters with real-time visualization:

- **Parameter panel** with tabs for body selection, simulation settings, material properties, and numerical options
- **Run manager** for tracking and comparing multiple simulation runs
- **Five plot types**: diurnal curves, depth profiles, heatmaps, flux, and combined views
- **Publication-quality export** (PDF, PNG, SVG) and CSV data export
- **Parameter sweeps** across any model variable (linear or log spacing)
- **Custom depth profiles** via the layer editor for non-standard regolith structures
- **JPL Horizons search** for finding any solar system body by name

![heat1d GUI with depth profile editor and diurnal temperature curves](docs/images/heat1d_screenshot.png)

Install with `pip install "heat1d[gui]"` and launch with `heat1d-gui`.

## C Implementation

The `c/` directory contains a standalone C implementation optimized for batch
runs and high performance. It supports all four solvers, reads the **same YAML
configuration files** as the Python version, and can model **any airless
planetary body** (Moon, Mercury, asteroids, icy satellites, etc.) by
configuring the physical and orbital parameters.

### Building

```bash
cd c/
brew install fftw libyaml   # macOS; see c/README.md for Linux packages
make
make test                    # 24 validation checks
```

### Running

**Legacy positional arguments** (backward compatible):

```bash
# Moon equator, TI=55, H=0.06, highland albedo
./heat1d 0 55 0.06 0.12 > temperature.txt

# Apollo 15 site (26 N, mare albedo, Fourier solver)
./heat1d 26 55 0.06 0.06 3 > apollo15.txt

# Implicit solver with adaptive timestepping
./heat1d 0 55 0.06 0.12 2 480 480 1.0 > implicit_adaptive.txt
```

**YAML configuration** (shares config files with Python):

```bash
# Use the Python example config
./heat1d --config ../python/heat1d/examples/moon_default.yaml --ti 55

# Override latitude and albedo for an Apollo site
./heat1d --config ../python/heat1d/examples/moon_default.yaml \
    --lat 26 --ti 55 --albedo 0.06

# Inspect the full configuration
./heat1d --config moon.yaml --verbose 2>config_dump.txt
```

Output is temperature data to stdout (one row per time step, one column per
depth layer), plus `loctime.txt` and `profile_z_dz_rho_k.txt` in the working
directory. See [c/README.md](c/README.md) for complete documentation.

## Theory

`heat1d` solves the 1-D heat equation in porous planetary regolith:

$$\rho \, c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial z} \left( K \frac{\partial T}{\partial z} \right)$$

with the following depth- and temperature-dependent material properties from Hayne et al. (2017):

- **Density**: exponential transition from surface (1100 kg/m³) to depth (1800 kg/m³) with scale height $H$
- **Thermal conductivity**: phonon (contact) conductivity plus a radiative $T^3$ component: $K = K_c \left[1 + \chi (T/350)^3 \right]$
- **Heat capacity**: 4th-order polynomial in temperature following Hemingway et al. (1981)
- **Surface boundary**: radiative equilibrium with angle-dependent albedo $A(\theta) = A_0 + a(\theta/45{^\circ})^3 + b(\theta/90{^\circ})^8$
- **Bottom boundary**: constant geothermal heat flux (0.018 W/m² for the Moon)

For full derivations, see the [Theory](python/docs/theory.md) documentation and [Hayne et al. (2017)](https://doi.org/10.1002/2017JE005387).

## Documentation

Detailed documentation is organized by topic:

| Topic | Description |
|---|---|
| [Theory](python/docs/theory.md) | Heat equation physics and flux formulation |
| [Thermophysical Properties](python/docs/properties.md) | Density, conductivity, and heat capacity models |
| [Boundary Conditions](python/docs/boundary.md) | Surface energy balance and bottom heat flux |
| [Numerical Methods](python/docs/numerical.md) | Solver equations, stability, and comparison |
| [Spatial Grid](python/docs/grid.md) | Non-uniform grid construction and skin depth |
| [Initialization](python/docs/initialization.md) | Temperature profile initialization |
| [Equilibration](python/docs/equilibration.md) | Convergence to periodic steady state |
| [Validation](python/docs/validation.md) | Comparison with Apollo and Diviner data |
| [CLI Reference](python/docs/cli.md) | Command-line options and examples |
| [API Reference](python/docs/api.rst) | Python class and function documentation |
| [C Implementation](c/README.md) | Building, running, and configuring the C code |
| [MATLAB Implementation](matlab/README.md) | MEX-file version for MATLAB |

Rendered documentation: [heat1d.readthedocs.io](https://heat1d.readthedocs.io)

## Validation

The model is validated against lunar temperature measurements from the Apollo 15 and 17 heat flow experiments and the Diviner Lunar Radiometer Experiment. All 8 validation checks pass within the published uncertainties across all four solvers.

```bash
heat1d --validate
```

See the [Validation](python/docs/validation.md) documentation for details.

## Citation

If you use `heat1d` in your research, please cite:

> Hayne, P. O., Bandfield, J. L., Siegler, M. A., Vasavada, A. R., Ghent, R. R., Williams, J.-P., Greenhagen, B. T., Aharonson, O., Elder, C. M., Lucey, P. G., & Paige, D. A. (2017). Global regolith thermophysical properties of the Moon from the Diviner Lunar Radiometer Experiment. *Journal of Geophysical Research: Planets*, 122, 2371-2400. [doi:10.1002/2017JE005387](https://doi.org/10.1002/2017JE005387)

## License

MIT License. See [LICENSE](LICENSE) for details.

## Authors

- **Paul O. Hayne** - University of Colorado Boulder (paul.hayne@lasp.colorado.edu)
- **K.-Michael Aye** - Package maintainer
