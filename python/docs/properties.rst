Thermophysical Properties
=========================

The regolith thermophysical properties in ``heat1d`` follow the models described in
Hayne et al. (2017), Eqs. A2--A6.

Density
-------

Bulk density increases exponentially with depth from a surface value :math:`\rho_s`
to a deep value :math:`\rho_d`:

.. math::
   \rho(z) = \rho_d - (\rho_d - \rho_s) \exp(-z / H)

where :math:`H` is the *H-parameter*, the e-folding scale depth.

For the Moon (Table A1 of Hayne et al., 2017):

=================  ====================  ==========
Parameter          Value                 Units
=================  ====================  ==========
:math:`\rho_s`     1100                  kg m\ :sup:`-3`
:math:`\rho_d`     1800                  kg m\ :sup:`-3`
:math:`H`          0.07                  m
=================  ====================  ==========

Contact Conductivity
--------------------

The *contact* (phonon) thermal conductivity follows the same depth profile:

.. math::
   K_c(z) = K_d - (K_d - K_s) \exp(-z / H)

=================  ===================  =========================
Parameter          Value                Units
=================  ===================  =========================
:math:`K_s`        7.4e-4               W m\ :sup:`-1` K\ :sup:`-1`
:math:`K_d`        3.4e-3               W m\ :sup:`-1` K\ :sup:`-1`
=================  ===================  =========================

Radiative Conductivity
----------------------

At elevated temperatures, radiative heat transfer between grains enhances the
effective thermal conductivity. The total conductivity is:

.. math::
   K = K_c (1 + \chi \, T^3 / 350^3)

where :math:`\chi` is a dimensionless parameter controlling the strength of
radiative conduction. At :math:`T = 350` K, the radiative contribution equals
:math:`\chi \cdot K_c`. The default value for the Moon is :math:`\chi = 2.7`
(Eq. A5 of Hayne et al., 2017).

Heat Capacity
-------------

The heat capacity :math:`c_p(T)` is a polynomial function of temperature, based on
laboratory data from Hemingway et al. (1981) and Ledlow et al. (1992):

.. math::
   c_p(T) = c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4

where the coefficients are stored in the ``planets`` package and are specific to
each planetary body. The polynomial yields non-physical (negative) values for
:math:`T < 1.3` K, but is valid for :math:`T \gtrsim 10` K.

Thermal Inertia
---------------

The thermal inertia is defined as:

.. math::
   I = \sqrt{K \rho c_p}

It controls the amplitude of diurnal temperature variations. Low thermal inertia
(loose regolith) produces large day-night contrasts, while high thermal inertia
(rock) produces small contrasts.

API Reference
-------------

.. automodule:: heat1d.properties
   :members:
   :undoc-members:
   :noindex:
