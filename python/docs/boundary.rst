Boundary Conditions
===================

The thermal model requires boundary conditions at the surface (top) and at depth
(bottom). These follow Hayne et al. (2017), Appendix A1.1 and A2.1
(Eqs. A7--A12, A21--A30).

Surface Boundary Condition
--------------------------

The surface energy balance determines the surface temperature :math:`T_s`:

.. math::
   \varepsilon \sigma T_s^4 = Q_s + K \left. \frac{\partial T}{\partial z} \right|_{z=0}

where:

- :math:`\varepsilon` is the infrared emissivity (0.95 for the Moon)
- :math:`\sigma` is the Stefan-Boltzmann constant
- :math:`Q_s` is the absorbed solar flux
- The last term is the conductive heat flux from below

The absorbed solar flux depends on the solar incidence angle :math:`\theta`:

.. math::
   Q_s = \frac{S_0}{r^2} (1 - A(\theta)) \cos\theta

where :math:`S_0` is the solar constant at 1 AU, :math:`r` is the heliocentric
distance in AU, and :math:`A(\theta)` is the incidence angle-dependent albedo.

Angle-Dependent Albedo
^^^^^^^^^^^^^^^^^^^^^^

The albedo model follows Keihm (1984) and Vasavada et al. (2012), as used in
Hayne et al. (2017), Eq. A10:

.. math::
   A(\theta) = A_0 + a \left(\frac{\theta}{\pi/4}\right)^3
   + b \left(\frac{\theta}{\pi/2}\right)^8

where :math:`A_0` is the normal bolometric Bond albedo (0.12 for highland,
0.07 for mare; Hayne et al. 2017), :math:`a = 0.06`, and :math:`b = 0.25`
for the Moon. The angle-dependent enhancement causes the effective hemispheric
Bond albedo to exceed :math:`A_0`.

Newton's Method
^^^^^^^^^^^^^^^

The surface energy balance is a nonlinear equation in :math:`T_s` and is solved
iteratively using Newton's root-finding method. The iteration:

.. math::
   T_s^{(k+1)} = T_s^{(k)} - \frac{f(T_s^{(k)})}{f'(T_s^{(k)})}

converges when :math:`|\Delta T| < \epsilon` (default: 0.1 K). A maximum
iteration count (default: 100) prevents infinite loops.

The surface temperature gradient is approximated using a second-order forward
difference:

.. math::
   \left. \frac{\partial T}{\partial z} \right|_{z=0} \approx
   \frac{-3 T_0 + 4 T_1 - T_2}{2 \Delta z_0}

Bottom Boundary Condition
-------------------------

The bottom boundary is set by the interior heat flux :math:`Q_b`:

.. math::
   T_N = T_{N-1} + \frac{Q_b}{K_{N-1}} \Delta z_{N-1}

For the Moon, :math:`Q_b = 0.018` W m\ :sup:`-2` (Langseth et al., 1976).

API Reference
-------------

.. automodule:: heat1d.boundary
   :members:
   :undoc-members:
   :noindex:
