# Numerical Methods

`heat1d` implements four numerical methods for the 1-D heat equation: three
finite-difference time-stepping schemes following Hayne et al. (2017), Appendix
A2 (Eqs. A13--A20), and a frequency-domain Fourier-matrix solver that eliminates
time-stepping entirely.

## Finite Difference Discretization

The heat equation on a non-uniform grid is discretized using the flux-conservative
form. For interior node $i$, the semi-discrete equation is:

$$
\rho_i c_{p,i} \frac{\partial T_i}{\partial t} =
p_i K_{i-1} (T_{i-1} - T_i) + q_i K_i (T_{i+1} - T_i)
$$

where the geometric coefficients $p_i$ and $q_i$ account for the
non-uniform grid spacing:

$$
p_i = \frac{2 \Delta z_i}{\Delta z_{i-1} \Delta z_i (\Delta z_{i-1} + \Delta z_i)}
$$

$$
q_i = \frac{2 \Delta z_{i-1}}{\Delta z_{i-1} \Delta z_i (\Delta z_{i-1} + \Delta z_i)}
$$

## Explicit Scheme (Forward Euler)

The explicit scheme (Eq. A17) updates temperatures using only values from the
current time step:

$$
T_i^{n+1} = T_i^n + \frac{\Delta t}{\rho_i c_{p,i}}
\left[ \alpha_i T_{i-1}^n - (\alpha_i + \beta_i) T_i^n + \beta_i T_{i+1}^n \right]
$$

where $\alpha_i = p_i K_{i-1}$ and $\beta_i = q_i K_i$.

**Stability**: The explicit scheme is conditionally stable, requiring the
CFL condition:

$$
\Delta t \leq F \cdot \min_i \frac{\rho_i c_{p,i} \Delta z_i^2}{K_i}
$$

where $F \leq 0.5$ is the Fourier mesh number. For typical Moon parameters,
this limits the time step to $\sim 3000$ s ($\sim 830$ steps per
lunar day).

## Implicit Scheme (Backward Euler)

The fully implicit scheme evaluates spatial derivatives at the new time level
$n+1$:

$$
-a_i T_{i-1}^{n+1} + (1 + a_i + b_i) T_i^{n+1} - b_i T_{i+1}^{n+1} = T_i^n
$$

where:

$$
a_i = \frac{\Delta t \, p_i K_{i-1}}{\rho_i c_{p,i}}, \quad
b_i = \frac{\Delta t \, q_i K_i}{\rho_i c_{p,i}}
$$

This forms a tridiagonal system that is solved using the Thomas algorithm
(see below).

**Stability**: The implicit scheme is *unconditionally stable*, meaning the time
step is limited by accuracy rather than stability. This allows much larger time
steps (e.g., $\sim 24$ steps per lunar day), providing a $\sim 35\times$
speedup over the explicit scheme.

**Accuracy**: First-order in time ($O(\Delta t)$), same as explicit.

## Crank-Nicolson Scheme (Semi-Implicit)

The Crank-Nicolson scheme averages the explicit and implicit contributions,
achieving second-order accuracy in time:

$$
-\frac{a_i}{2} T_{i-1}^{n+1} + \left(1 + \frac{a_i + b_i}{2}\right) T_i^{n+1}
- \frac{b_i}{2} T_{i+1}^{n+1} =
\frac{a_i}{2} T_{i-1}^n + \left(1 - \frac{a_i + b_i}{2}\right) T_i^n
+ \frac{b_i}{2} T_{i+1}^n
$$

The left-hand side forms a tridiagonal system with half-coefficients, while the
right-hand side uses the old temperature values.

**Stability**: Unconditionally stable, like the fully implicit scheme.

**Accuracy**: Second-order in time ($O(\Delta t^2)$), the most accurate
of the three schemes for a given time step.

**Boundary values**: The Crank-Nicolson scheme requires both old and new
boundary values. The old values (before the BC update) contribute to the explicit
half of the RHS, while the new values contribute to the implicit half of the LHS.

## Thomas Algorithm (TDMA)

Both the implicit and Crank-Nicolson schemes produce tridiagonal linear systems
of the form:

$$
\mathbf{A} \mathbf{x} = \mathbf{d}
$$

where $\mathbf{A}$ is tridiagonal with sub-diagonal $\ell$, main
diagonal $m$, and super-diagonal $u$.

The Thomas algorithm (Tri-Diagonal Matrix Algorithm) solves this in $O(n)$
time using forward elimination and back-substitution:

**Forward sweep**:

$$
c'_0 = \frac{u_0}{m_0}, \quad
d'_0 = \frac{d_0}{m_0}
$$

$$
c'_i = \frac{u_i}{m_i - \ell_{i-1} c'_{i-1}}, \quad
d'_i = \frac{d_i - \ell_{i-1} d'_{i-1}}{m_i - \ell_{i-1} c'_{i-1}}
$$

**Back-substitution**:

$$
x_n = d'_n, \quad
x_i = d'_i - c'_i x_{i+1}
$$

The algorithm is numerically stable for the heat equation because the coefficient
matrix is always diagonally dominant: $|m_i| = 1 + a_i + b_i > |a_i| + |b_i|$.

## Fourier-Matrix Solver (Frequency Domain)

The Fourier-matrix solver takes a fundamentally different approach: instead of
marching forward in time, it solves for the periodic steady-state temperature
directly in the frequency domain. This eliminates time-stepping and equilibration
entirely, providing ~100--1000× speedup over the finite-difference methods.

### Principle

For a periodic surface forcing with period $P$, the temperature at any depth
can be decomposed into Fourier harmonics:

$$
T(z, t) = \bar{T}(z) + \sum_{n=1}^{N/2} \hat{T}_n(z) \, e^{i n \omega_0 t} + \text{c.c.}
$$

where $\omega_0 = 2\pi / P$ is the fundamental angular frequency and
$\bar{T}(z)$ is the time-mean (DC) temperature profile. Each harmonic
propagates independently through the subsurface, with amplitude decaying and
phase shifting according to the thermal properties of each layer.

### Transmission Matrices

For a single frequency $\omega$, the temperature and heat flux at the top and
bottom of a homogeneous layer of thickness $d$ are related by a 2×2 transmission
matrix (analogous to electrical transmission lines):

$$
\begin{pmatrix} \hat{T} \\ \hat{q} \end{pmatrix}_{\text{top}} =
\begin{pmatrix} \cosh(qd) & \frac{\sinh(qd)}{kq} \\
kq \sinh(qd) & \cosh(qd) \end{pmatrix}
\begin{pmatrix} \hat{T} \\ \hat{q} \end{pmatrix}_{\text{bottom}}
$$

where $q = \sqrt{i\omega / \kappa}$ is the complex thermal wavenumber, $k$ is
the thermal conductivity, and $\kappa = k / (\rho c_p)$ is the thermal
diffusivity. The matrix product over all layers gives the global transfer from
the surface to the bottom boundary, and the **surface thermal impedance** is:

$$
Z_{\text{surf}}(\omega) = \frac{\hat{T}_{\text{surf}}}{\hat{q}_{\text{surf}}} = \frac{P_{00}}{P_{10}}
$$

where $P_{00}$ and $P_{10}$ are elements of the cumulative matrix product.

### Nonlinear Surface Radiation

The surface energy balance $\varepsilon \sigma T_s^4 = Q_s + q_{\text{cond}}$
is nonlinear in $T_s$. The solver handles this via Newton iteration in the
time domain. The conductive heat flux at the surface is computed from the
frequency-domain admittance (inverse impedance) as a circulant matrix $\mathbf{C}$:

$$
\mathbf{q}_{\text{cond}} = \mathbf{C} \, \mathbf{T}_s
$$

where $\mathbf{C}$ is constructed from the inverse FFT of the admittance
spectrum $Y(\omega_n) = 1 / Z_{\text{surf}}(\omega_n)$. The Newton update at
each iteration is:

$$
\delta \mathbf{T}_s = \left( \mathbf{J} + \mathbf{C} \right)^{-1} \mathbf{R}
$$

where $\mathbf{J} = \text{diag}(4 \varepsilon \sigma T_s^3)$ is the Jacobian of
the radiation term and $\mathbf{R}$ is the residual. This converges in 5--15
iterations.

### Thermal Pumping (Rectification)

The nonlinear temperature dependence of thermal conductivity
($K \propto 1 + \chi T^3/350^3$) produces a net downward heat transport known
as *thermal pumping* or the *solid-state greenhouse effect*. Because conductivity
is higher when the near-surface is hot (daytime), more heat flows downward during
the day than upward at night. This elevates subsurface temperatures above what
a linear model would predict.

The solver captures this through an outer iteration loop:

1. **Freeze properties** at the current equilibrium profile $\bar{T}(z)$
2. **Solve** the linearized frequency-domain problem (inner Newton loop)
3. **Compute rectification flux** $J_{\text{pump}}(z) = \langle k(T) \, \partial T'/\partial z \rangle$ from the time-domain reconstruction of $T(z,t)$ and the exact nonlinear $k(T)$
4. **Update the equilibrium profile** by integrating $d\bar{T}/dz = (Q_b - J_{\text{pump}}) / \langle k \rangle$ downward from the surface
5. **Repeat** until the mean surface temperature converges (typically 3--5 outer iterations)

### Depth Reconstruction

Once the surface temperature is converged, the full $T(z, t)$ field is
reconstructed using the depth transfer functions:

$$
\hat{T}_n(z) = \hat{T}_{n,\text{surf}} \cdot \frac{P_{00}(z)}{P_{00,\text{surf}}}
$$

The DC component uses the equilibrium profile $\bar{T}(z)$ (which includes the
rectification correction), and the AC components are propagated from the surface
via the transmission matrices. An inverse FFT recovers $T(z, t)$.

### Performance

The Fourier-matrix solver is ~1000× faster than time-stepping for a single
diurnal cycle because:

- No equilibration orbits are required (the solution is already periodic)
- All frequencies are solved simultaneously via matrix operations
- The Newton iteration converges in $O(10)$ iterations rather than $O(10^3)$ time steps
- A complete lunar diurnal cycle is solved in ~100 ms

### Limitations

- Assumes periodic forcing (cannot model transient events like eclipses)
- Linearizes subsurface conduction around the equilibrium profile (corrected by
  outer iteration, but may be less accurate for extreme temperature swings)
- Does not support PSR (bowl-shaped crater) geometry

## Solver Comparison

| Scheme | Method | Accuracy | Stability | Steps/lunar day | Relative Speed |
|---|---|---|---|---|---|
| Explicit | Forward Euler | $O(\Delta t)$ | Conditional | ~830 | 1× |
| Implicit | Backward Euler + TDMA | $O(\Delta t)$ | Unconditional | ~24 | ~35× |
| Crank-Nicolson | Semi-implicit + TDMA | $O(\Delta t^2)$ | Unconditional | ~24 | ~35× |
| Fourier-matrix | Frequency domain | Spectral | N/A (periodic) | N/A | ~1000× |

For time-stepping applications, the Crank-Nicolson scheme is recommended: it
offers second-order accuracy with the same unconditional stability as the fully
implicit scheme.

For periodic steady-state problems (the most common use case), the **Fourier-matrix
solver** is strongly preferred. It is the default equilibration solver
(`equil_solver = "fourier-matrix"` in `Configurator`), and when used as the
primary solver (`solver = "fourier-matrix"`), it bypasses equilibration entirely.

### Fourier-Matrix as Equilibration Solver

Even when the output phase uses a time-stepping solver (e.g., for eclipse
modeling or external flux series), the Fourier-matrix solver is used by default
for the equilibration phase. It computes the periodic steady-state temperature
profile directly, then initializes the time-stepping solver from $T(t=0, z)$
(local noon). This eliminates the need for multi-orbit spin-up and ensures the
time-stepping solver starts from a well-converged state. See
[Equilibration](equilibration.md) for details.
