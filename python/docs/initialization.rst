Initialization
==============

Temperature initialization follows Hayne et al. (2017), Appendix A2.3
(Eqs. A34--A36).

Initial Temperature Profile
----------------------------

The temperature profile is initialized to the *equilibrium mean temperature*,
which approximates the expected time-averaged temperature for a rapidly rotating
body:

.. math::
   T_{eq}(\phi) = \frac{T_{rad}(\phi)}{\sqrt{2}}

where :math:`T_{rad}` is the radiative equilibrium temperature at local noon:

.. math::
   T_{rad}(\phi) = \left[ \frac{(1 - A) S_0 \cos\phi}{\varepsilon \sigma} \right]^{1/4}

and :math:`\phi` is the latitude.

This initialization provides a reasonable starting point that speeds up
convergence to the periodic steady state. All layers are initially set to the
same temperature.

Property Initialization
-----------------------

After the temperature profile is set, the heat capacity and thermal conductivity
profiles are computed from the initial temperatures:

1. Heat capacity: :math:`c_p(T)` via the polynomial fit
2. Thermal conductivity: :math:`K(z, T)` combining depth-dependent contact
   conductivity and temperature-dependent radiative conductivity
