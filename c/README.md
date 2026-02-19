# heat1d -- C Implementation

Native C implementation of the 1-D thermal model for planetary regolith.
This is the most computationally efficient version, solving the heat equation
with depth- and temperature-dependent thermophysical properties following
[Hayne et al. (2017)](https://doi.org/10.1002/2017JE005387).

Although the default parameters are calibrated for the Moon, all physical and
orbital properties are configurable via YAML or CLI arguments, making the code
suitable for **any airless planetary body** (Mercury, asteroids, icy satellites,
etc.).

## Prerequisites

External dependencies:
- [FFTW3](http://www.fftw.org/) -- used by the Fourier-matrix solver
- [libyaml](https://pyyaml.org/wiki/LibYAML) -- used by the YAML configuration parser

LAPACK/BLAS are also required on Linux; macOS uses the built-in Accelerate framework.

**macOS** (Homebrew):

```bash
brew install fftw libyaml
```

**Ubuntu / Debian**:

```bash
sudo apt install libfftw3-dev liblapack-dev libyaml-dev
```

**Fedora / RHEL**:

```bash
sudo dnf install fftw-devel lapack-devel libyaml-devel
```

## Building

From the `c/` directory:

```bash
make
```

This produces the `heat1d` executable. The Makefile auto-detects macOS vs.
Linux for library paths and framework linkage.

To build and run the validation test suite:

```bash
make test
```

## Quick Start

Run the Moon at the equator with default parameters:

```bash
./heat1d 0 55 0.06 0.12 > temperature.txt
```

This simulates one diurnal cycle at the lunar equator with:
- Latitude = 0 degrees
- Thermal inertia = 55 (standard regolith at 273 K)
- H-parameter = 0.06 m (density/conductivity scale depth)
- Albedo = 0.12 (highland)

Temperature data (one row per time step, one column per depth layer) is written
to stdout, so redirect it to a file. Two auxiliary files are written to the
current directory:
- `loctime.txt` -- local solar time (hours past noon) for each output row
- `profile_z_dz_rho_k.txt` -- depth grid: z [m], dz [m], density [kg/m³], conductivity [W/m/K]

### More Examples

```bash
# Apollo 15 site (26°N, mare albedo, Fourier solver)
./heat1d 26 55 0.06 0.06 3 > apollo15.txt

# High latitude with Crank-Nicolson solver
./heat1d 60 55 0.06 0.12 1 > lat60.txt

# Low thermal inertia (fluffy regolith)
./heat1d 0 30 0.06 0.12 > fluffy.txt
```

## YAML Configuration

The C code can read the **same YAML configuration files** used by the Python
implementation. This allows all physical and numerical parameters to be set
without recompiling.

```bash
./heat1d --config ../python/heat1d/examples/moon_default.yaml --ti 55
```

The `--config` flag loads the YAML file, then optional CLI flags override
individual values:

| Flag | Description |
|---|---|
| `--config <file>` | YAML configuration file (required for YAML mode) |
| `--lat <degrees>` | Override latitude |
| `--ti <value>` | Override thermal inertia |
| `--H <value>` | Override H-parameter |
| `--albedo <value>` | Override albedo |
| `--flux <file>` | External flux file |
| `--verbose` | Print full configuration to stderr |

### YAML File Format

The YAML file has three main sections: `planet`, `numerical`, and `physical`.
All fields are optional; unspecified values use compile-time defaults.

```yaml
planet:
  name: "Moon"
  albedo: 0.12
  emissivity: 0.95
  ks: 7.4e-4          # surface conductivity [W/m/K]
  kd: 3.4e-3          # deep conductivity [W/m/K]
  rhos: 1100.0         # surface density [kg/m^3]
  rhod: 1800.0         # deep density [kg/m^3]
  H: 0.06              # H-parameter [m]
  Qb: 0.018            # interior heat flow [W/m^2]
  day: 2550240.0       # synodic solar day [s]
  obliquity: 0.0       # axial tilt [rad]
  eccentricity: 0.0
  rAU: 1.0             # semi-major axis [AU]

latitude: 0.0            # degrees
ndays: 1                 # number of output days
solver: "implicit"       # "explicit", "crank-nicolson", "implicit", or "fourier-matrix"

numerical:
  accuracy: 1.0                # adaptive timestepping tolerance [K]
  output_interval: 53130       # output spacing [s]
  fourier_number: 0.5
  layers_per_skin_depth: 10
  layer_growth_factor: 5
  skin_depths_to_bottom: 20
  surface_temp_accuracy: 0.1   # [K]
  equilibration_years: 1

physical:
  solar_constant: 1361.0       # [W/m^2]
  radiative_conductivity: 2.7  # chi parameter
```

## Command-Line Arguments (Legacy Mode)

```
./heat1d <lat> <T.I.> <H> <albedo> [solver] [equil_nperday] [nperday_output] [adaptive_tol] [flux_file]
```

The first four arguments are required; the rest are optional.

| # | Argument | Units | Description | Default / Typical |
|---|---|---|---|---|
| 1 | `lat` | degrees | Latitude | -- (required) |
| 2 | `T.I.` | J m⁻² K⁻¹ s⁻¹/² | Thermal inertia at 273 K | 55 (standard regolith) |
| 3 | `H` | m | Density/conductivity e-folding scale depth | 0.06 |
| 4 | `albedo` | -- | Normal bolometric Bond albedo | 0.12 (highland), 0.06 (mare) |
| 5 | `solver` | -- | 0 = explicit, 1 = Crank-Nicolson, 2 = implicit, 3 = Fourier | 0 |
| 6 | `equil_nperday` | -- | Time steps per day during equilibration | 480 |
| 7 | `nperday_output` | -- | Output samples per day | 480 |
| 8 | `adaptive_tol` | K | Adaptive step-doubling tolerance (0 = off) | 0 |
| 9 | `flux_file` | -- | Path to external flux file (overrides solar flux) | none |

### Solvers

| Code | Solver | Notes |
|---|---|---|
| 0 | **Explicit** (Forward Euler) | CFL-limited, ~830 steps/day. Simple and robust. |
| 1 | **Crank-Nicolson** | Unconditionally stable, 2nd-order in time. Recommended for time-stepping. |
| 2 | **Implicit** (Backward Euler) | Unconditionally stable, 1st-order in time. |
| 3 | **Fourier-matrix** | Frequency-domain solver. Solves periodic steady state directly (~1000x faster). Requires FFTW3. |

## Output Format

**stdout** (redirect to a file):

Each row is one output time step. Each column is one depth layer (column 0 =
surface). Values are temperature in Kelvin, space-separated:

```
388.47 385.12 378.93 ... 251.84
387.21 384.98 378.81 ... 251.84
...
```

**loctime.txt**:

One value per line, the local solar time in hours past noon (0 = noon, 6 =
sunset, 12 = midnight, 18 = sunrise, 24 = noon again):

```
0.00
0.05
0.10
...
```

**profile_z_dz_rho_k.txt**:

One row per layer: depth [m], layer thickness [m], density [kg/m³], phonon
conductivity [W/m/K]:

```
0.0000 4.256e-04 1100 7.400e-04
0.0004 5.107e-04 1104 7.432e-04
...
```

## Plotting Results

### With Python

```python
import numpy as np
import matplotlib.pyplot as plt

T = np.loadtxt("temperature.txt")
lt = np.loadtxt("loctime.txt")

plt.plot(lt, T[:, 0])
plt.xlabel("Local Time (hours past noon)")
plt.ylabel("Surface Temperature [K]")
plt.title("Lunar Equatorial Surface Temperature")
plt.show()
```

### With gnuplot

```gnuplot
set xlabel "Local Time (hours past noon)"
set ylabel "Surface Temperature [K]"
plot "< paste loctime.txt temperature.txt" using 1:2 with lines title "Surface"
```

## Physical Parameters

All physical constants are defined in `heat1dfun.h`. Key values for the Moon
(Table A1 of Hayne et al., 2017):

| Parameter | Symbol | Value | Units |
|---|---|---|---|
| Solar constant | S₀ | 1361 | W m⁻² |
| Stefan-Boltzmann | σ | 5.67 × 10⁻⁸ | W m⁻² K⁻⁴ |
| IR emissivity | ε | 0.95 | -- |
| Synodic period | P | 2.55024 × 10⁶ | s |
| Surface density | ρₛ | 1100 | kg m⁻³ |
| Deep density | ρ_d | 1800 | kg m⁻³ |
| Surface conductivity | Kₛ | 7.4 × 10⁻⁴ | W m⁻¹ K⁻¹ |
| Deep conductivity | K_d | 3.4 × 10⁻³ | W m⁻¹ K⁻¹ |
| Radiative parameter | χ | 2.7 | -- |
| Geothermal heat flux | Q | 0.018 | W m⁻² |

## Source Files

| File | Description |
|---|---|
| `heat1d.c` | Main program: CLI/YAML parsing, grid setup, simulation driver |
| `heat1dfun.c` | Core solver: time-stepping, boundary conditions, material properties |
| `heat1dfun.h` | Constants, structures (`layerT`, `profileT`), function prototypes |
| `orbitfun.c` | Orbital mechanics: solar zenith angle, hour angle, orbit update |
| `orbitfun.h` | Orbital function prototypes |
| `fourier_solver.c` | Frequency-domain Fourier-matrix solver (requires FFTW3) |
| `fourier_solver.h` | Fourier solver structures and prototypes |
| `yaml_config.c` | YAML configuration parser (requires libyaml) |
| `yaml_config.h` | Configuration struct (`configT`) and prototypes |
| `test_validate.c` | Validation tests against Hayne et al. (2017) Table A2 |
| `Makefile` | Build system (auto-detects macOS/Linux) |

## Reference

Hayne, P. O., et al. (2017). Global regolith thermophysical properties of the
Moon from the Diviner Lunar Radiometer Experiment. *Journal of Geophysical
Research: Planets*, 122, 2371-2400.
[doi:10.1002/2017JE005387](https://doi.org/10.1002/2017JE005387)
