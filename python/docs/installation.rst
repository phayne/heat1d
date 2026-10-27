.. highlight:: shell

============
Installation
============

.. note::

   The ``heat1d`` package on PyPI is an obsolete release (0.3.2) and is no
   longer updated. Do not use ``pip install heat1d``; install from GitHub
   as described below.

heat1d requires Python 3.9 or newer.


From GitHub (recommended)
-------------------------

The Python package lives in the ``python/`` subdirectory of the `Github repo`_.
Install the latest version directly with pip:

.. code-block:: console

    $ pip install "heat1d @ git+https://github.com/phayne/heat1d.git#subdirectory=python"

To include the graphical user interface (PySide6):

.. code-block:: console

    $ pip install "heat1d[gui] @ git+https://github.com/phayne/heat1d.git#subdirectory=python"

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.


From sources
------------

For development, clone the repository and do an editable install:

.. code-block:: console

    $ git clone https://github.com/phayne/heat1d.git
    $ cd heat1d/python
    $ pip install -e ".[dev,gui]"

Alternatively, download the `tarball`_, extract it, and run
``pip install .`` from its ``python/`` subdirectory.


.. _pip: https://pip.pypa.io
.. _Python installation guide: https://docs.python-guide.org/starting/installation/
.. _Github repo: https://github.com/phayne/heat1d
.. _tarball: https://github.com/phayne/heat1d/tarball/master
