# Equilibration

Before producing output, the thermal model must reach a periodic steady state
where the temperature profile repeats from one diurnal cycle to the next.

## Approach

The model runs for a user-specified number of equilibration orbits
(`NYEARSEQ`, default: 1) before recording output. During equilibration, the
temperature profile evolves from the initial guess toward the periodic steady
state.

## Convergence

The depth to which temperatures equilibrate scales as:

$$
z_{eq} \sim z_s \sqrt{N_{cycles}}
$$

where $z_s$ is the skin depth and $N_{cycles}$ is the number of
cycles completed. Deep subsurface temperatures require more cycles to converge.

For the Moon:

- Surface temperatures equilibrate within 1--2 lunar days
- Subsurface temperatures at 1 m depth may require 5+ orbits for convergence
  to < 1 K accuracy

## Computational Cost

The equilibration phase dominates the total computation time. The solver choice
has a large impact:

| Scheme | Steps per day | Relative cost |
|---|---|---|
| Explicit | ~830 | 1.0× (baseline) |
| Crank-Nicolson | ~24 | ~0.03× |
| Implicit | ~24 | ~0.03× |

The implicit and Crank-Nicolson schemes achieve $\sim 35\times$ speedup
because they are not constrained by the CFL stability limit. For long
equilibration runs (many orbits or annual cycles), the implicit solvers are
strongly recommended.

## Configuration

Equilibration parameters are set in the `Configurator` class or YAML config:

- `NYEARSEQ`: Number of equilibration orbits (default: 1)
- `equil_dt`: Equilibration timestep in seconds (default: `None` = `day/48`)
- `DTBOT`: Bottom temperature convergence criterion in K (default: 0.1)
