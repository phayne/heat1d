# Theory

The 1-D heat equation governs the evolution of temperature $T$ as a function
of depth $z$ and time $t$:

$$
\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}
\left( K \frac{\partial T}{\partial z} \right)
$$

where $\rho$ is the bulk density, $c_p$ is the specific heat capacity,
and $K$ is the thermal conductivity.

For planetary regolith, all three material properties depend on depth, and both
$c_p$ and $K$ also depend on temperature, making the equation
nonlinear.

## Physical Interpretation

The heat equation expresses conservation of energy. The left-hand side is the
rate of thermal energy storage per unit volume. The right-hand side is the
divergence of the conductive heat flux $q = -K \, \partial T / \partial z$.

In a planetary context, the upper boundary is driven by the diurnal cycle of
solar insolation and thermal emission to space, while the lower boundary is
governed by interior heat flow. The model propagates surface temperature
variations downward, with the depth of penetration controlled by the *thermal
skin depth* (see [Spatial Grid](grid.md)).

## Flux Form

For numerical solution, the heat equation is written in flux-conservative form.
The heat flux across layer boundaries is:

$$
q_{i+1/2} = -K_{i+1/2} \frac{T_{i+1} - T_i}{\Delta z_{i+1/2}}
$$

where $K_{i+1/2}$ is evaluated at the interface between layers $i$
and $i+1$. This form ensures energy conservation on the discrete grid
(see [Numerical Methods](numerical.md) for discretization details).

## Reference

Hayne, P. O., et al. (2017). Appendix A1: "Thermal model".
*J. Geophys. Res. Planets*, 122, 2371--2400.
