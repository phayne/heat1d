import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

mpl.rcParams["font.size"] = 18
mpl.rcParams["lines.linewidth"] = 2
mpl.rcParams["lines.linestyle"] = "--"
mpl.rcParams["figure.constrained_layout.use"] = True


def profile_plot(model, ax1=None):
    m = model
    if ax1 is None:
        fig, ax1 = plt.subplots()
    else:
        fig = ax1.figure

    # Depth profiles of min, mean, max temperature
    # ax1.set_xlim(-10+np.nanmin(m.T[:,0]),10+np.nanmax(m.T[:,0]))
    ax1.set_yscale("log")
    ax1.plot(m.T.max(0), m.profile.z, label="$T_\mathrm{max}$")
    ax1.plot(m.T.min(0), m.profile.z, label="$T_\mathrm{min}$")
    ax1.plot(m.T.mean(0), m.profile.z, ls="-", label="$T_\mathrm{avg}$")
    ax1.set_xlabel("Temperature, $T$ (K)")
    ax1.set_ylabel("Depth, $z$ (m)")
    ax1.legend(frameon=False)
    # ax1.set_ylim(np.nanmax(m.profile.z), m.profile.z[1])
    ax1.set_ylim(1.5, m.profile.z[1])


def diurnal_curves(model, ax2=None):
    m = model
    if ax2 is None:
        fig, ax2 = plt.subplots()
    else:
        fig = ax2.figure

    # Diurnal temperature curves
    ax2.set_xlim(m.lt.min(), np.nanmax(m.lt) + 40)
    ax2.set_ylim(-5 + m.T.min(), m.T.max() + 5)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    for i, z in enumerate(m.profile.z):
        ax2.plot(
            m.lt,
            m.T[:, i],
            ls="-",
            lw=1,
            label="{:.3f}".format(z),
            color=plt.cm.magma(i * 10),
        )
    ax2.legend(frameon=False, title="Depth (m):", fontsize=8)
    ax2.set_xlabel("Local Time (hr past noon)")
    ax2.set_ylabel("Temperature, $T$ (K)")


def plot_profile_and_diurnals(model, save=False):
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(14, 5))
    profile_plot(model, ax1)
    diurnal_curves(model, ax2)
    fig.suptitle(f"Chi: {model.profile.chi}")
    if save:
        fig.savefig("heat1d_example.pdf")
