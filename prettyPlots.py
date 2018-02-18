"""
=============
Pretty Plots
=============

This module provides methods for making prettier plots,
with custom annotation

KNOWN ISSUE: Incompatible w/ pandas

Author: Paul O. Hayne

"""

import matplotlib
import numpy as np

rcParams = matplotlib.rcParams
rc = matplotlib.rc

# Set overall style of plots:
# Thicker lines, larger text, LaTeX rendering
def setStyle():
    setFont(16)
    rc('lines', linewidth=2, markersize=8)
    rc('axes', titlesize=22, labelsize=20)
    rc('xtick', labelsize=18)
    rc('ytick', labelsize=18)

# Font settings
def setFont(font_size):
    font = {'family' : 'serif',
            'size'   : font_size}
    rc('font', **font)
    rc('text', usetex=True)

#####################
# degreeLabelFormat #
#####################
# Format text with degree symbol
@np.vectorize
def degreeLabelFormat(x):
    if rcParams['text.usetex'] and not rcParams['text.latex.unicode']:
        return r"$%0.0f^\circ$" % x
    else:
        return "%0.0f\u00b0" % x