======
heat1d
======


.. image:: https://readthedocs.org/projects/heat1d/badge/?version=latest
        :target: https://heat1d.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Thermal model for planetary science applications (Python version)


* Free software: MIT license
* Documentation: https://heat1d.readthedocs.io.

See the notebook_ for an example.

.. _notebook: https://nbviewer.jupyter.org/github/phayne/heat1d/blob/master/python/notebooks/heat1d_example.ipynb


Features
--------

* Modeling 1D heat transfer in planetary surfaces
* Provides library of required planetary surface parameters
* Plotting helpers (see notebook_)

Installation
------------

The ``heat1d`` release on PyPI is obsolete and no longer updated. Install the
current version from GitHub (requires Python 3.9 or newer):

.. code-block:: bash

   pip install "heat1d @ git+https://github.com/phayne/heat1d.git#subdirectory=python"

Add the ``[gui]`` extra (``"heat1d[gui] @ git+..."``) for the graphical user
interface. See the `installation docs`_ for development installs.

.. _`installation docs`: https://heat1d.readthedocs.io/en/latest/installation.html


Credits
---------

Author: P. Hayne

Package maintainer: K.-Michael Aye

This package was created with Cookiecutter_ and the forked `michaelaye/cookiecutter-pypackage-conda`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`michaelaye/cookiecutter-pypackage-conda`: https://github.com/michaelaye/cookiecutter-pypackage-conda
