Spatial Grid
============

The finite difference grid uses geometrically increasing layer thicknesses,
following Hayne et al. (2017), Appendix A2.2 (Eqs. A31--A33).

Thermal Skin Depth
------------------

The *thermal skin depth* :math:`z_s` is the characteristic depth of penetration
of a periodic temperature wave:

.. math::
   z_s = \sqrt{\frac{\kappa P}{\pi}}

where :math:`P` is the forcing period (e.g., the diurnal period) and
:math:`\kappa = K / (\rho c_p)` is the thermal diffusivity.

For the Moon with surface properties, :math:`z_s \approx 4\text{--}5` cm.

Grid Construction
-----------------

The grid starts at :math:`z = 0` (the surface) with an initial layer thickness:

.. math::
   \Delta z_0 = \frac{z_s}{m}

where :math:`m` is the number of layers within the first skin depth (default: 10).

Layer thickness grows geometrically with depth:

.. math::
   \Delta z_{i+1} = \Delta z_i \left(1 + \frac{1}{n}\right)

where :math:`n` controls the growth rate (default: 5). Larger :math:`n` gives
more uniform layers; smaller :math:`n` gives faster growth.

The grid extends to a total depth of :math:`b` skin depths (default: 20), ensuring
that the bottom boundary is far enough below the surface that diurnal temperature
variations are negligible.

Grid Parameters
---------------

===========  ========  =========================================
Parameter    Default   Description
===========  ========  =========================================
:math:`m`    10        Layers per skin depth
:math:`n`    5         Growth factor (dz[i+1] = dz[i]*(1+1/n))
:math:`b`    20        Total depth in skin depths
===========  ========  =========================================

For the Moon, these defaults produce approximately 45 layers extending to a
depth of about 1 m.

API Reference
-------------

.. automodule:: heat1d.grid
   :members:
   :undoc-members:
   :noindex:
