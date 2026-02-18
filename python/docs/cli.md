# Command-Line Interface

`heat1d` provides a command-line interface for running thermal models with
optional YAML configuration files.

## Basic Usage

```bash
# Default: Moon equator, explicit solver, 1 day
heat1d

# Specify latitude and solver
heat1d --lat 45 --solver implicit

# Use a YAML config file
heat1d moon.yaml

# Override YAML values with CLI flags
heat1d moon.yaml --lat 30 --solver crank-nicolson

# Run validation suite
heat1d --validate

# Quiet mode (no progress output)
heat1d --quiet --no-plot
```

## Options

```text
Usage: heat1d [OPTIONS] [CONFIG_FILE]

Arguments:
  CONFIG_FILE    YAML configuration file (optional)

Options:
  --lat FLOAT          Latitude [degrees]
  --ndays INTEGER      Number of output days
  --solver CHOICE      explicit | crank-nicolson | implicit
  --planet TEXT         Planet name (e.g., Moon)
  --chi FLOAT          Radiative conductivity parameter
  --albedo FLOAT       Surface albedo
  --output-dir PATH    Output directory
  --prefix TEXT        Output filename prefix
  --no-plot            Suppress plot generation
  --validate           Run Moon validation suite
  --quiet              Suppress progress output
  -h, --help           Show help and exit
```

## YAML Configuration

The CLI accepts YAML configuration files that specify all model parameters.
Command-line options override YAML values when both are provided.

Example YAML (`moon_default.yaml`):

```yaml
planet:
  name: "Moon"
  # Override any planet property:
  # albedo: 0.12
  # ks: 7.4e-4

latitude: 0.0
ndays: 1

solver: "explicit"

numerical:
  fourier_number: 0.5
  layers_per_skin_depth: 10
  layer_growth_factor: 5
  skin_depths_to_bottom: 20
  steps_per_day: 24
  surface_temp_accuracy: 0.1
  bottom_temp_accuracy: 0.1
  equilibration_years: 1

physical:
  solar_constant: 1361.0
  radiative_conductivity: 2.7

output:
  directory: "output/"
  prefix: "moon_eq"
```

### YAML Sections

**planet**: Planet name and optional property overrides. Properties that can be
overridden: `albedo`, `emissivity`, `ks`, `kd`, `rhos`, `rhod`,
`H`, `cp0`, `Qb`.

**numerical**: Grid and solver parameters mapped to `Configurator` fields.

**physical**: Physical constants (solar constant, radiative conductivity).

**output**: Output directory and filename prefix.

## Output Files

The CLI generates the following output files:

- `{prefix}_temperature.csv`: Temperature vs. local time for all depths
- `{prefix}_grid.csv`: Depth grid with density and conductivity
- `{prefix}_plot.png`: Combined profile and diurnal curve plot (unless `--no-plot`)
