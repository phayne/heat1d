# Theory

`heat1d` solves the 1-D heat equation for planetary regolith with depth- and
temperature-dependent thermophysical properties, following the formulation of
Hayne et al. (2017), Appendix A (Eqs. A1--A36).

## Governing Equation

The 1-D heat equation governs the evolution of temperature $T$ as a function
of depth $z$ and time $t$:

$$
\rho \, c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial z} \left( K \frac{\partial T}{\partial z} \right) \qquad \text{(A1)}
$$

where $\rho(z)$ is the bulk density, $c_p(T)$ is the specific heat capacity,
and $K(z, T)$ is the thermal conductivity. All three properties depend on depth,
and $c_p$ and $K$ also depend on temperature, making the equation nonlinear.

The left-hand side is the rate of thermal energy storage per unit volume. The
right-hand side is the divergence of the conductive heat flux
$q = -K \, \partial T / \partial z$. The model propagates surface temperature
variations downward, with the depth of penetration controlled by the *thermal
skin depth* $z_s = \sqrt{\kappa P / \pi}$, where $\kappa = K / (\rho c_p)$ is
the thermal diffusivity and $P$ is the forcing period
(see [Spatial Grid](grid.md)).

## Thermophysical Properties

The regolith properties follow depth and temperature profiles calibrated to
lunar data (Hayne et al. 2017, Table A1). For full details, see
[Thermophysical Properties](properties.md).

### Density

Bulk density increases exponentially from a surface value $\rho_s$ to a
deep value $\rho_d$ (Eq. A2):

$$
\rho(z) = \rho_d - (\rho_d - \rho_s) \, e^{-z/H} \qquad \text{(A2)}
$$

where $H$ is the *H-parameter*, the e-folding scale depth. For the Moon,
$\rho_s = 1100$ kg m⁻³, $\rho_d = 1800$ kg m⁻³, and $H \approx 0.06$ m.

### Thermal Conductivity

The total thermal conductivity combines a phonon (contact) term $K_c$ and a
radiative $T^3$ term (Eqs. A3--A5):

$$
K = K_c \left[ 1 + \chi \left(\frac{T}{350}\right)^3 \right] \qquad \text{(A4)}
$$

The contact conductivity follows the same depth profile as density:

$$
K_c(z) = K_d - (K_d - K_s) \frac{\rho_d - \rho(z)}{\rho_d - \rho_s} \qquad \text{(A5)}
$$

For the Moon, $K_s = 7.4 \times 10^{-4}$ W m⁻¹ K⁻¹,
$K_d = 3.4 \times 10^{-3}$ W m⁻¹ K⁻¹, and $\chi = 2.7$.

The temperature dependence has an important physical consequence: because
conductivity is higher when the near-surface is hot (daytime), more heat flows
downward during the day than upward at night. This net downward *thermal
pumping* — the **solid-state greenhouse effect** — elevates subsurface
temperatures above what a linear model would predict.

### Heat Capacity

The specific heat capacity is a polynomial in temperature (Eq. A6):

$$
c_p(T) = c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4 \qquad \text{(A6)}
$$

based on data from Hemingway et al. (1981) and Ledlow et al. (1992). Valid for
$T \gtrsim 10$ K. See [Thermophysical Properties](properties.md) for
coefficient values.

### Thermal Inertia

The thermal inertia $I = \sqrt{K \rho c_p}$ controls the amplitude of diurnal
temperature variations. Low $I$ (loose, porous regolith) produces large
day-night temperature contrasts; high $I$ (compacted soil or rock) produces
small contrasts. Typical lunar regolith has
$I \approx 55\ \mathrm{J\,m^{-2}\,K^{-1}\,s^{-1/2}}$ at 273 K.

## Surface Boundary Condition

The surface temperature $T_s$ is determined by the energy balance between
absorbed solar radiation, thermal emission, and conduction into the subsurface
(Eq. A7):

$$
\varepsilon \, \sigma \, T_s^4 = Q_s + K \left. \frac{\partial T}{\partial z} \right|_{z=0} \qquad \text{(A7)}
$$

where $\varepsilon$ is the infrared emissivity ($0.95$ for the Moon), $\sigma$
is the Stefan-Boltzmann constant, and $Q_s$ is the absorbed solar flux. This
equation states that the outgoing thermal radiation (left side) must equal the
sum of absorbed sunlight and heat conducted from below (right side).

At night, $Q_s = 0$ and surface cooling is balanced solely by subsurface
conduction. The surface cools rapidly at first, then asymptotically as the
conducted flux diminishes, producing the characteristic plateau shape of the
nighttime cooling curve.

### Absorbed Solar Flux

The absorbed flux on a horizontal surface at latitude $\phi$ is (Eq. A9):

$$
Q_s(t) = \frac{S_0}{r^2} \left(1 - A(\theta)\right) \cos\theta \qquad \text{(A9)}
$$

where $S_0 = 1361$ W m⁻² is the solar constant at 1 AU, $r$ is the
heliocentric distance in AU, $A(\theta)$ is the albedo (which depends on
incidence angle $\theta$), and $\cos\theta$ gives the projected area factor.

The solar incidence angle $\theta$ for a horizontal surface depends on
latitude $\phi$, solar declination $\delta$, and hour angle $h = 2\pi t / P$:

$$
\cos\theta = \sin\phi \sin\delta + \cos\phi \cos\delta \cos h \qquad \text{(A11)}
$$

The hour angle $h$ is measured from local noon, and the flux is clipped to zero
when $\cos\theta < 0$ (the Sun is below the horizon):

$$
\psi(x) = \frac{1}{2}\left(\cos x + \lvert\cos x\rvert\right) \qquad \text{(A10)}
$$

For a body with orbital eccentricity $e$ and obliquity $\epsilon$, the solar
declination angle $\delta$ varies over the orbital period according to the
relationship $\sin\delta = \sin\epsilon \sin L_s$, where $L_s$ is the
areocentric longitude.

### Angle-Dependent Albedo

The Bond albedo increases with incidence angle following Keihm (1984) and
Vasavada et al. (2012):

$$
A(\theta) = A_0 + a \left(\frac{\theta}{\pi/4}\right)^3 + b \left(\frac{\theta}{\pi/2}\right)^8 \qquad \text{(A8)}
$$

where $A_0$ is the normal-incidence albedo (0.12 for highland, 0.06 for mare),
$a = 0.06$, and $b = 0.25$ for the Moon. The cubic and octic terms cause a
rapid increase in albedo at grazing incidence angles, enhancing reflection near
sunrise and sunset. The effective hemispheric Bond albedo (integrated over all
angles) is significantly higher than $A_0$ — approximately 0.23 for
$A_0 = 0.12$.

### Newton's Method for Surface Temperature

Because Eq. A7 is nonlinear in $T_s$ (through the $T_s^4$ emission term and the
$T^3$-dependent conductivity), it is solved iteratively using Newton's method
(Eqs. A21--A29).

The function whose root is sought is:

$$
f(T_s) = \varepsilon \, \sigma \, T_s^4 - Q_s - K_0 \frac{-3T_0 + 4T_1 - T_2}{2 \Delta z_0} \qquad \text{(A25)}
$$

where the surface temperature gradient uses a second-order forward difference
(Eq. A24). Its derivative with respect to $T_s$ is:

$$
f'(T_s) = 4 \varepsilon \sigma T_s^3 - 3 B_0 T_s^2 \frac{4T_1 - 3T_0 - T_2}{2 \Delta z_0} + \frac{3}{2 \Delta z_0} \left(K_{c,0} + B_0 T_s^3\right) \qquad \text{(A29)}
$$

where $B_0$ is the radiative conductivity prefactor at the surface. The
iteration $T_s^{(k+1)}$ = $T_s^{(k)} - f/f'$ converges when
$|\Delta T| < \epsilon$ (default: 0.1 K), typically in 2--5 iterations.

## Bottom Boundary Condition

The lower boundary applies a constant geothermal heat flux $Q_b$ (Eq. A12):

$$
\left. \frac{\partial T}{\partial z} \right|_{z=z_*} = \frac{Q_b}{K} \qquad \text{(A12)}
$$

In the finite-difference discretization, this becomes (Eq. A30):

$$
T_N = T_{N-1} + \frac{Q_b}{K_{N-1}} \Delta z_{N-1} \qquad \text{(A30)}
$$

For the Moon, $Q_b = 0.018$ W m⁻² (Langseth et al., 1976), corresponding
to a temperature gradient of a few K m⁻¹ at depth. This is a small but
non-negligible correction: the heat flux maintains a weak upward temperature
gradient that slightly elevates subsurface temperatures.

The bottom boundary is placed at a depth of $\sim 20$ thermal skin depths below
the surface ($\sim 1$ m for the Moon). At this depth, diurnal temperature
variations are attenuated by a factor of $e^{-20} \approx 2 \times 10^{-9}$,
so the fixed-flux condition does not influence the diurnal surface temperature.
See [Spatial Grid](grid.md) for details on grid construction.

## Flux Form for Numerical Solution

For numerical solution, the heat equation is written in flux-conservative form.
The heat flux across layer boundaries is (Eq. A15):

$$
q_{i+1/2} \approx K_i \frac{T_{i+1} - T_i}{\Delta z_i} \qquad \text{(A15)}
$$

The flux gradient at node $i$ uses the fluxes on either side (Eq. A16):

$$
\left.\frac{\partial q}{\partial z}\right|_i \approx \frac{2}{\Delta z_i + \Delta z_{i-1}} \left[ K_i \frac{T_{i+1} - T_i}{\Delta z_i} - K_{i-1} \frac{T_i - T_{i-1}}{\Delta z_{i-1}} \right] \qquad \text{(A16)}
$$

This form ensures energy conservation on the non-uniform grid. The geometric
coefficients $p_i$ and $q_i$ (Eq. A18) absorb the grid-spacing factors:

$$
p_i = \frac{2 \Delta z_i}{\Delta z_{i-1} \Delta z_i (\Delta z_{i-1} + \Delta z_i)}, \quad q_i = \frac{2 \Delta z_{i-1}}{\Delta z_{i-1} \Delta z_i (\Delta z_{i-1} + \Delta z_i)}
$$

so the temperature update simplifies to:

$$
T_i^{n+1} = T_i^n + \frac{\Delta t}{\rho_i c_{p,i}} \left[ \alpha_i T_{i-1}^n - (\alpha_i + \beta_i) T_i^n + \beta_i T_{i+1}^n \right] \qquad \text{(A17)}
$$

where

$$
\alpha_i = p_i K_{i-1}, \quad \beta_i = q_i K_i
$$

For the full discretization of all four solvers (explicit, implicit,
Crank-Nicolson, Fourier-matrix), see [Numerical Methods](numerical.md).

## Temperature Initialization

To reduce equilibration time, the initial temperature profile is set using
analytic estimates (Eqs. A34--A36).

**Surface**: radiative equilibrium at local noon:

$$
T_0 = \left(\frac{(1 - A) \, F_\odot}{\varepsilon \, \sigma}\right)^{1/4} \qquad \text{(A34)}
$$

**Bottom**: isothermal body equilibrium:

$$
T_N = \frac{T_0}{\sqrt{2}} \qquad \text{(A35)}
$$

**Subsurface**: exponential interpolation:

$$
T_i = T_N - (T_N - T_0) \, e^{-z_i/H} \qquad \text{(A36)}
$$

This gives a plausible starting profile. The model then equilibrates to the
self-consistent periodic steady state over several orbital cycles. When using
the Fourier-matrix solver for equilibration (the default), convergence is
achieved directly without time-stepping. See [Equilibration](equilibration.md)
and [Initialization](initialization.md) for details.

## Summary of Key Parameters (Moon)

| Parameter | Symbol | Value | Units | Source |
|---|---|---|---|---|
| Solar constant | $S_0$ | 1361 | W m⁻² | Kopp and Lean (2011) |
| Stefan-Boltzmann | $\sigma$ | $5.67 \times 10^{-8}$ | W m⁻² K⁻⁴ | |
| IR emissivity | $\varepsilon$ | 0.95 | -- | Logan et al. (1972) |
| Synodic period | $P$ | $2.55 \times 10^{6}$ | s | |
| Surface density | $\rho_s$ | 1100 | kg m⁻³ | Hayne et al. (2013) |
| Deep density | $\rho_d$ | 1800 | kg m⁻³ | Carrier et al. (1991) |
| Surface conductivity | $K_s$ | $7.4 \times 10^{-4}$ | W m⁻¹ K⁻¹ | Hayne et al. (2017) |
| Deep conductivity | $K_d$ | $3.4 \times 10^{-3}$ | W m⁻¹ K⁻¹ | Hayne et al. (2017) |
| Radiative parameter | $\chi$ | 2.7 | -- | Vasavada et al. (2012) |
| H-parameter | $H$ | 0.06 | m | Hayne et al. (2017) |
| Normal albedo (highland) | $A_0$ | 0.12 | -- | Vasavada et al. (2012) |
| Albedo coefficient | $a$ | 0.06 | -- | Hayne et al. (2017) |
| Albedo coefficient | $b$ | 0.25 | -- | Hayne et al. (2017) |
| Geothermal heat flux | $Q_b$ | 0.018 | W m⁻² | Langseth et al. (1976) |

## Reference

Hayne, P. O., Bandfield, J. L., Siegler, M. A., Vasavada, A. R., Ghent, R. R.,
Williams, J.-P., Greenhagen, B. T., Aharonson, O., Elder, C. M., Lucey, P. G.,
& Paige, D. A. (2017). Global regolith thermophysical properties of the Moon
from the Diviner Lunar Radiometer Experiment. *J. Geophys. Res. Planets*, 122,
2371--2400. [doi:10.1002/2017JE005387](https://doi.org/10.1002/2017JE005387)
