Theory
======

The 1-D heat equation governs the evolution of temperature :math:`T` as a function
of depth :math:`z` and time :math:`t`:

.. math::
   \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}
   \left( K \frac{\partial T}{\partial z} \right)

where :math:`\rho` is the bulk density, :math:`c_p` is the specific heat capacity,
and :math:`K` is the thermal conductivity.

For planetary regolith, all three material properties depend on depth, and both
:math:`c_p` and :math:`K` also depend on temperature, making the equation
nonlinear.

Physical Interpretation
-----------------------

The heat equation expresses conservation of energy. The left-hand side is the
rate of thermal energy storage per unit volume. The right-hand side is the
divergence of the conductive heat flux :math:`q = -K \, \partial T / \partial z`.

In a planetary context, the upper boundary is driven by the diurnal cycle of
solar insolation and thermal emission to space, while the lower boundary is
governed by interior heat flow. The model propagates surface temperature
variations downward, with the depth of penetration controlled by the *thermal
skin depth* (see :doc:`grid`).

Flux Form
---------

For numerical solution, the heat equation is written in flux-conservative form.
The heat flux across layer boundaries is:

.. math::
   q_{i+1/2} = -K_{i+1/2} \frac{T_{i+1} - T_i}{\Delta z_{i+1/2}}

where :math:`K_{i+1/2}` is evaluated at the interface between layers :math:`i`
and :math:`i+1`. This form ensures energy conservation on the discrete grid
(see :doc:`numerical` for discretization details).

Reference
---------

Hayne, P. O., et al. (2017). Appendix A1: "Thermal model".
*J. Geophys. Res. Planets*, 122, 2371--2400.
