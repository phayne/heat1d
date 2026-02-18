Validation
==========

The ``heat1d`` model is validated against lunar temperature data from
Hayne et al. (2017), Table A2, using constraints from Diviner observations
and Apollo heat flow experiments.

Validation Constraints
----------------------

The following constraints are used (Table A2 of Hayne et al., 2017):

**Equatorial temperatures** (latitude = 0, highland albedo A\ :sub:`h` = 0.12):

=============================  ===========  ===========
Constraint                     Value (K)    Tolerance
=============================  ===========  ===========
Peak noon temperature          385          +/- 5
Midnight temperature           101          +/- 5
Minimum nighttime temperature  95           +/- 5
=============================  ===========  ===========

**Apollo site mean temperatures** (mare albedo A\ :sub:`h` = 0.06):

====================================  ==========  ==========  ===========
Constraint                            Latitude    Value (K)   Tolerance
====================================  ==========  ==========  ===========
Apollo 15 surface mean T              26 N        211         +/- 5
Apollo 15 subsurface (0.83 m) mean T  26 N        252         +/- 5
Apollo 17 surface mean T              20 N        216         +/- 5
Apollo 17 subsurface (0.13 m) mean T  20 N        256         +/- 5
====================================  ==========  ==========  ===========

Mare vs. Highland Albedo
~~~~~~~~~~~~~~~~~~~~~~~~

The equator checks use the default Moon highland normal bolometric Bond albedo
(A\ :sub:`h` = 0.12). The Apollo landing sites are in dark mare regions with
significantly lower albedo. Following Hayne et al. (2017), which reports
A\ :sub:`h` = 0.12 for highland and A\ :sub:`h` = 0.07 for mare, the Apollo
checks use a mare albedo of A\ :sub:`h` = 0.06 appropriate for the particularly
dark basaltic floors at Hadley Rille (Apollo 15) and Taurus-Littrow (Apollo 17).

Density/Conductivity Scale Depth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The validation suite uses the Hayne et al. (2017) Table A1 standard value for
the density and conductivity e-folding scale depth (*H* = 0.06 m), overriding
the ``planets`` package default of *H* = 0.07 m. This improves the fit to
Diviner nighttime cooling observations at all latitudes (RMS residual drops
from ~1 K to ~0.3 K).

Energy Conservation
-------------------

In addition to temperature comparisons, the validation suite checks energy
conservation by computing the stored energy change over one diurnal cycle.
For a well-equilibrated model, the relative energy imbalance should be < 1%.

Running the Validation Suite
----------------------------

From the command line::

    heat1d --validate

Or from Python::

    from heat1d.validation import run_validation_suite
    results = run_validation_suite(solver="explicit", nyearseq=5)

The validation suite generates four plots:

1. **Diurnal equator curve** -- Surface temperature vs. local time at the equator,
   with reference values marked
2. **Multi-latitude diurnal curves** -- Surface temperature at 0, 15, 30, 45, 60,
   and 75 degrees latitude
3. **Nighttime cooling curves** -- Surface temperature during the lunar night at
   multiple latitudes with Diviner regolith temperature observations, similar
   to Figure A2
4. **Mean temperature vs. latitude** -- Diurnal mean surface and subsurface
   temperature vs. latitude, with Apollo 15 and 17 reference data and error bars
   for both surface and subsurface measurements

Validation Results
------------------

With the Hayne et al. (2017) Table A1 standard properties (highland albedo for
equator, mare albedo for Apollo sites, *H* = 0.06 m), all 8 validation checks
pass. The Apollo checks use a finer grid (m=20, b=30) and longer equilibration
(25 orbits) to ensure the 0.83 m subsurface temperature is well converged:

.. code-block:: text

    [PASS] equator_peak_noon_T: 388.5 K (ref: 385.0 +/- 5.0 K)
    [PASS] equator_midnight_T: 100.2 K (ref: 101.0 +/- 5.0 K)
    [PASS] equator_min_night_T: 93.7 K (ref: 95.0 +/- 5.0 K)
    [PASS] energy_conservation: relative error = 0.0000
    [PASS] apollo15_surface_mean_T: 209.1 K (ref: 211.0 +/- 5.0 K)
    [PASS] apollo15_subsurface_mean_T: 252.9 K (ref: 252.0 +/- 5.0 K)
    [PASS] apollo17_surface_mean_T: 211.7 K (ref: 216.0 +/- 5.0 K)
    [PASS] apollo17_subsurface_mean_T: 255.4 K (ref: 256.0 +/- 5.0 K)
    8/8 checks passed

API Reference
-------------

.. automodule:: heat1d.validation
   :members:
   :undoc-members:
   :noindex:
