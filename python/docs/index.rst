heat1d Documentation
====================

**heat1d** is a one-dimensional thermal model for planetary science applications.
It solves the heat equation in a layered, porous regolith using finite differences,
with temperature-dependent thermal properties and a non-uniform spatial grid.

This documentation follows the structure of the Appendix from:

    Hayne, P. O., et al. (2017). Global regolith thermophysical properties of the
    Moon from the Diviner Lunar Radiometer Experiment. *Journal of Geophysical
    Research: Planets*, 122, 2371--2400.
    `doi:10.1002/2017JE005387 <https://doi.org/10.1002/2017JE005387>`_

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   readme
   installation
   usage
   cli

.. toctree::
   :maxdepth: 2
   :caption: Theory & Methods

   theory
   properties
   boundary
   numerical
   grid
   initialization
   equilibration

.. toctree::
   :maxdepth: 2
   :caption: Validation

   validation

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   authors
   history

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
