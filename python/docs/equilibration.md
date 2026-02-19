# Equilibration

Before producing output, the thermal model must reach a periodic steady state
where the temperature profile repeats from one diurnal cycle to the next.
`heat1d` provides two approaches: **Fourier-matrix equilibration** (default,
instantaneous) and **time-stepping equilibration** (iterative, for special cases).

## Fourier-Matrix Equilibration (Default)

By default, equilibration uses the Fourier-matrix solver
(`equil_solver = "fourier-matrix"` in `Configurator`). This computes the
periodic steady-state temperature $T(z, t)$ directly in the frequency domain,
without any time-stepping or multi-orbit convergence. The Fourier solution is
exact for the linearized problem (plus outer iterations for nonlinear
conductivity; see [Numerical Methods](numerical.md)).

The process:

1. The Fourier-matrix solver computes the complete diurnal cycle $T(z, t)$
2. The profile is initialized from $T(t=0, z)$ (local noon)
3. Material properties ($K$, $c_p$) are updated to match the initialized temperatures
4. The time-stepping solver begins output from this converged state

This eliminates the traditional multi-orbit spin-up entirely and is the
recommended approach for all standard simulations. The equilibration phase
completes in ~100 ms regardless of how many output days are requested.

**When Fourier equilibration is used**: always, unless `equil_solver` is
explicitly changed. This includes all time-stepping output solvers (explicit,
implicit, Crank-Nicolson). When `solver = "fourier-matrix"`, equilibration is
not needed at all because the solver itself produces the periodic steady state.

## Time-Stepping Equilibration

For specialized applications (e.g., validating the Fourier solver against
time-stepping, or modeling situations where the Fourier solver's periodicity
assumption is inappropriate), the equilibration solver can be set to a
time-stepping scheme:

```python
config = Configurator(solver="explicit", equil_solver="implicit")
```

In this mode, the model runs for `NYEARSEQ` orbits before recording output.
During equilibration, the temperature profile evolves from the initial guess
toward the periodic steady state through repeated diurnal cycles.

### Convergence

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

### Computational Cost

When using time-stepping equilibration, the equilibration phase dominates the
total computation time. The solver choice has a large impact:

| Scheme | Steps per day | Relative cost |
|---|---|---|
| Fourier-matrix (default) | N/A | ~0.001× |
| Explicit | ~830 | 1.0× (baseline) |
| Crank-Nicolson | ~24 | ~0.03× |
| Implicit | ~24 | ~0.03× |

The Fourier-matrix solver achieves ~1000× speedup over explicit time-stepping
because it computes the periodic steady state directly without iterating through
individual diurnal cycles. For time-stepping equilibration, the implicit and
Crank-Nicolson schemes achieve ~35× speedup because they are not
constrained by the CFL stability limit.

## Configuration

Equilibration parameters are set in the `Configurator` class or YAML config:

- `equil_solver`: Solver for equilibration phase (default: `"fourier-matrix"`)
- `NYEARSEQ`: Number of equilibration orbits for time-stepping (default: 1)
- `equil_dt`: Equilibration timestep in seconds (default: `None` = `day/48`)
- `DTBOT`: Bottom temperature convergence criterion in K (default: 0.1)
