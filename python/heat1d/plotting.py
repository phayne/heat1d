"""
Plotting module

This module provides the old style settings from the former prettyPlots but makes use of the newer matplotlib stylesheets for increased flexibility.

Authors: Paul O. Hayne (setStyle, setFont), Michael Aye (the rest)
"""
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import numpy as np


# look at plt.style.available for other styles
plt.style.use("seaborn")

# mpl.rcParams["figure.constrained_layout.use"] = True
# plt.rcParams.update(
#     {
#         "text.usetex": True,
#         "font.family": "serif",
#         "font.serif": ["Palatino"],
#         "lines.linewidth": 2,
#         "lines.markersize": 8,
#     }
# )


def profile_axes_labels(ax):
    ax.set_xlabel("Temperature, $T$ [K]")
    ax.set_ylabel("Depth, $z$ [m]")


def profile_plot(model, ax=None):
    m = model
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Depth profiles of min, mean, max temperature
    # ax1.set_xlim(-10+np.nanmin(m.T[:,0]),10+np.nanmax(m.T[:,0]))
    ax.set_yscale("log")
    ax.plot(m.T.max(0), m.profile.z, label="$T_\mathrm{max}$")
    ax.plot(m.T.min(0), m.profile.z, label="$T_\mathrm{min}$")
    ax.plot(m.T.mean(0), m.profile.z, ls="-", label="$T_\mathrm{avg}$")
    profile_axes_labels(ax)
    ax.legend(frameon=False)
    # ax1.set_ylim(np.nanmax(m.profile.z), m.profile.z[1])
    ax.set_ylim(1.5, m.profile.z[1])
    ax.set_title(f"Chi: {model.profile.chi}")


def diurnal_curves(model, ax=None):
    m = model
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Diurnal temperature curves
    ax.set_xlim(m.lt.min(), np.nanmax(m.lt) + 40)
    ax.set_ylim(-5 + m.T.min(), m.T.max() + 5)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    for i, z in enumerate(m.profile.z):
        ax.plot(
            m.lt,
            m.T[:, i],
            ls="-",
            lw=1,
            label="{:.3f}".format(z),
            color=plt.cm.magma(i * 10),
        )
    ax.legend(frameon=False, title="Depth (m):", fontsize=8)
    ax.set_xlabel("Local Time (hr past noon)")
    ax.set_ylabel("Temperature, $T$ [K]")
    ax.set_title(f"Chi: {model.profile.chi}")


def plot_last_surface_cooling(model, ax=None):
    m = model
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    m.last_surface_cooling.plot(ax=ax)
    ax.set_xlabel("Local Time (hr past noon)")
    ax.set_ylabel("Temperature, $T$ [K]")
    ax.set_title(f"Chi: {model.profile.chi}")


def plot_profile_and_diurnals(model, save=False):
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(9, 5))
    profile_plot(model, ax1)
    diurnal_curves(model, ax2)
    if save:
        fig.savefig("heat1d_example.pdf")


def compare_model_profiles(m1, m2):
    fig, ax = plt.subplots()
    ax.plot(m1.T.mean(0), m1.profile.z, ls="-", label=f"chi: {m1.profile.chi}")
    plt.plot(m2.T.mean(0), m2.profile.z, ls="-", label=f"chi: {m2.profile.chi}")
    ax.legend()
    profile_axes_labels(ax)
    ax.set_yscale("log")
    ax.set_ylim(1.5, m1.profile.z[1])
    ax.set_title("Average T-profile for different chi values")


def setStyle():
    """
    Set overall style of plots:
    Thicker lines, larger text, LaTeX rendering
    """
    import matplotlib

    rcParams = matplotlib.rcParams
    rc = matplotlib.rc
    setFont(rc, 16)
    rc("lines", linewidth=2, markersize=8)
    rc("axes", titlesize=22, labelsize=20)
    rc("xtick", labelsize=18)
    rc("ytick", labelsize=18)


def setFont(rc, font_size):
    font = {"family": "serif", "size": font_size}
    rc("font", **font)
    rc("text", usetex=True)


@np.vectorize
def degreeLabelFormat(x):
    "Format text with degree symbol."
    if rcParams["text.usetex"] and not rcParams["text.latex.unicode"]:
        return r"$%0.0f^\circ$" % x
    else:
        return "%0.0f\u00b0" % x
