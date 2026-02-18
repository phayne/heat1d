#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# heat1d documentation build configuration file

import os
import sys

# Insert project root so Sphinx can find the heat1d package
sys.path.insert(0, os.path.abspath(".."))

import heat1d

# -- General configuration ------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]

source_suffix = {
    ".rst": "restructuredtext",
}

master_doc = "index"

project = "heat1d"
copyright = "2024, Paul O. Hayne"
author = "Paul O. Hayne"

version = heat1d.__version__
release = heat1d.__version__

exclude_patterns = ["_build", "requirements.txt"]

pygments_style = "sphinx"

# -- Options for HTML output ----------------------------------------------

try:
    import sphinx_rtd_theme
    html_theme = "sphinx_rtd_theme"
except ImportError:
    html_theme = "alabaster"

html_static_path = ["_static"]

htmlhelp_basename = "heat1ddoc"

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {}

latex_documents = [
    ("index", "heat1d.tex", "heat1d Documentation", "Paul O. Hayne", "manual"),
]

man_pages = [
    ("index", "heat1d", "heat1d Documentation", ["Paul O. Hayne"], 1)
]

# -- Options for autodoc --------------------------------------------------

autodoc_member_order = "bysource"
