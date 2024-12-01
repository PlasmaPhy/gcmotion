r"""
Plots the time evolution of the dynamical variables theta, zeta, psi ,psip,
rho, Ptheta, Pzeta with respect to time.
"""

import numpy as np
import matplotlib.pyplot as plt
from re import split, IGNORECASE

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import ParticleEvolution as config
from gcmotion.entities.particle import Particle


def particle_evolution(particle: Particle, **args):
    r"""Plots the time evolution of all the dynamical variables and
    canonical momenta from the particle's calculated orbit.

    :meta public:

    Parameters
    ----------
    particle: :py:class:`~gcmotion.entities.particle.Particle`
        The current working particle.
    which : str, optional
        Which of the dynamical variables to show. Must be a string containing
        the variable names "theta", "zeta", "psi", "psip", "rho", "Ptheta",
        "Pzeta", seperated by "," or spaces. For example: "theta zeta psi,
        psip, Ptheta Pzeta". Defaults to "all".
    percentage : int, optional
        The percentage of the orbit to be plotted. Defaults to 100.
    units : str, optional
        The unit system. Can be either 'SI' or 'NU'. Defauls to "SI".
    """
    # Unpack parameters
    which = args.get("which", config.which)
    units = args.get("units", config.units)
    percentage = args.get("percentage", config.percentage)

    # This suffix is used to yank the correct attributes from particle.
    suffix = "NU" if units == "NU" else "" if units == "SI" else ""
    logger.info(
        "Plotting time evolutions in "
        + f"{"NU" if suffix == "NU" else "SI"}..."
    )

    # Make sure percentage is a valid number
    if percentage < 1 or percentage > 100:
        percentage = 100
        logger.warning("Invalid percentage: Plotting the whole thing...")
    points = int(np.floor(particle.t_eval.shape[0] * percentage / 100) - 1)

    # Parse "which" and yank the respective attributes form particle
    which = format_which(which)
    # The returned dict looks something like this:
    # which = ["theta", "zeta", "rho", "Pzeta"]
    # Now only the variables we want to plot appear, and in the correct order

    # Yank the trimmed time arrays from particle and store them in a dict
    t = getattr(particle, "t_eval" + suffix)[:points]
    variables = {}
    for key in which:
        if key in ["theta", "zeta"]:
            variables[key] = getattr(particle, key)[:points]
        else:
            variables[key] = getattr(particle, key + suffix)[:points]

    # Create figure
    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
        "facecolor": config.facecolor,
    }
    fig = plt.figure(**fig_kw)
    fig.suptitle(
        "Time evolution of dynamical variables",
        size=config.titlesize,
        color=config.titlecolor,
    )

    # Create requested axis (axeses?)
    # stack a subplot for every remaining variable
    ax = fig.subplots(len(which), 1, sharex=True)

    # Iterate over variables and axes and draw
    # Also define "Pzeta_idx" which stores the ax index of the Pzeta plot, we
    # need it later
    scatter_kw = {
        "s": config.s,
        "color": config.color,
        "marker": config.marker,
    }

    i = 0
    Pzeta_idx = None
    for key, value in variables.items():
        if key == "Pzeta":
            Pzeta_idx = i
        ax[i].scatter(t, value, **scatter_kw)
        i += 1

    # Label setup
    # iterate over axes and change each label
    # The label strings are set up at the end of this module
    labels = get_labels(particle, suffix)
    i = 0
    for key, value in variables.items():
        ax[i].set_ylabel(
            ylabel=labels[key],
            size=config.labelsize,
            labelpad=config.labelpad,
            rotation=90,
        )
        i += 1

    # Zoom out Pzeta y-axis
    if "Pzeta" in variables.keys():
        current_ylim = np.array(ax[Pzeta_idx].get_ylim())
        ax[Pzeta_idx].set_ylim(
            np.sort([current_ylim[0] / 3, current_ylim[1] * 3])
        )

    # Move yticks to the right
    for i in range(len(fig.axes)):
        ax[i].yaxis.tick_right()

    # Add common xlabel
    fig.supxlabel(f"t [{t.units:~P}]", size=config.labelsize)

    # Time axis limit
    ax[-1].set_xlim([t[0], t[-1]])
    plt.show()

    logger.info("--> Time evolutions successfully plotted.")


def format_which(which: str):
    r"""Returns a dictionary with the values to be plotted."""

    if which == "all":
        which = "theta zeta psi psip rho Ptheta Pzeta"

    # split spaces
    # Resulting list contains all the "words" delimited by " ", "," and "." and
    # enpyt strings, so we have to remove them.
    which_list = split(r" |,", which, flags=IGNORECASE)
    filtered = list(filter(bool, which_list))

    # Make sure invalid values are removed by make
    # This also makes sure it is in the correct order.
    res = {}
    res["theta"] = True if "theta" in filtered else False
    res["zeta"] = True if "zeta" in filtered else False
    res["psi"] = True if "psi" in filtered else False
    res["psip"] = True if "psip" in filtered else False
    res["rho"] = True if "rho" in filtered else False
    res["Ptheta"] = True if "Ptheta" in filtered else False
    res["Pzeta"] = True if "Pzeta" in filtered else False

    # Now res should be something like
    # res = {
    #     "theta": True,
    #     "zeta": True,
    #     "psi": True,
    #     "psip": False,
    #     "rho": True,
    #     "Ptheta": True,
    #     "Pzeta": True,
    # }

    # remove all the key-value pairs with value=False:
    result = {}
    for key, value in res.items():
        if value:
            result[key] = value

    return list(result.keys())


def get_labels(particle: Particle, suffix: str):
    r"""Sets up the ylabels since they are a bit bulky."""

    return {
        "theta": r"$\theta$",
        "zeta": r"$\zeta$",
        "psi": (
            r"$\psi$" + "\n" + f"[{getattr(particle, "psi"+suffix).units:~P}]"
        ),
        "psip": (
            r"$\psi_p$"
            + "\n"
            + f"[{getattr(particle, "psip"+suffix).units:~P}]"
        ),
        "rho": (
            r"$\rho$" + "\n" + f"[{getattr(particle, "rho"+suffix).units:~P}]"
        ),
        "Ptheta": (
            r"$P_\theta$"
            + "\n"
            + f"[{getattr(particle, "Ptheta"+suffix).units:~P}]"
        ),
        "Pzeta": (
            r"$P_\zeta$"
            + "\n"
            + f"[{getattr(particle, "Pzeta"+suffix).units:~P}]"
        ),
    }
