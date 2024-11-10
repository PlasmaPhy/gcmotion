r"""
Plots the time evolution of the dynamical variables
:math:`\theta, \zeta, \psi, \psi_p, \rho_{||}, P_\theta, P_\zeta`
with respect to time.

The time axis can be either [NU] (normalized units) of time, or 
seconds.
The percentage of the shown time evolution can be any integer 
from 1 to 100.

Example
-------

.. code-block:: python

    gcm.time_evolution(cwp, percentage=20, units="SI")

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import time_evolution as config


def time_evolution(cwp, percentage: int = 100, units: str = "SI"):
    r"""Plots the time evolution of all the dynamical variables and
    canonical momenta.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    percentage : int, optional
        The percentage of the orbit to be plotted. Defaults to 100.
    units : str, optional
        The unit system. Can be either 'SI' or 'NU'. Defauls to "SI".
    """
    suffix = "NU" if units == "NU" else "" if units == "SI" else ""
    logger.info(f"Plotting time evolutions in {"NU" if suffix=="NU" else "SI"}...")  # fmt: skip

    # Get all needed attributes first
    theta = getattr(cwp, "theta").copy()
    zeta = getattr(cwp, "zeta").copy()
    t = getattr(cwp, "t_eval" + suffix).copy()
    psi = getattr(cwp, "psi" + suffix).copy()
    rho = getattr(cwp, "rho" + suffix).copy()
    psip = getattr(cwp, "psip" + suffix).copy()
    Ptheta = getattr(cwp, "Ptheta" + suffix).copy()
    Pzeta = getattr(cwp, "Pzeta" + suffix).copy()

    if percentage < 1 or percentage > 100:
        percentage = 100
        print("Invalid percentage. Plotting the whole thing.")
        logger.warning("Invalid percentage: Plotting the whole thing...")

    points = int(np.floor(t.shape[0] * percentage / 100) - 1)

    # Plotting
    fig, ax = plt.subplots(7, 1, **config["fig_parameters"])
    ax[0].set_title("Time evolution of dynamical variables", c="b")
    ax[5].set_title("Time evolution of canonical momenta", c="b")

    # fmt: off
    ax[0].scatter(t[:points], theta[:points],   **config["scatter_args"])
    ax[1].scatter(t[:points], zeta[:points],    **config["scatter_args"])
    ax[2].scatter(t[:points], psi[:points],     **config["scatter_args"])
    ax[3].scatter(t[:points], psip[:points],    **config["scatter_args"])
    ax[4].scatter(t[:points], rho[:points],     **config["scatter_args"])
    ax[5].scatter(t[:points], Ptheta[:points],  **config["scatter_args"])
    ax[6].scatter(t[:points], Pzeta[:points],   **config["scatter_args"])

    lp = config["labelpad"]
    lp = lp*2 if suffix == "NU" else lp # give them NUs a bit more space
    loc = config["loc"]
    ax[0].set_ylabel(r"$\theta(t)$"   +"\n"+ rf"[${theta.units:P~}$]",  loc=loc, labelpad=lp, **config["ylabel_args"])
    ax[1].set_ylabel(r"$\zeta(t)$"    +"\n"+ rf"[${zeta.units:P~}$]",   loc=loc, labelpad=lp, **config["ylabel_args"])
    ax[2].set_ylabel(r"$\psi(t)$"     +"\n"+ rf"[${psi.units:P~}$]",    loc=loc, labelpad=lp, **config["ylabel_args"])
    ax[3].set_ylabel(r"$\psi_p(t)$"   +"\n"+ rf"[${psip.units:P~}$]",   loc=loc, labelpad=lp, **config["ylabel_args"])
    ax[4].set_ylabel(r"$\rho_{||}(t)$"+"\n"+ rf"[${rho.units:P~}$]",    loc=loc, labelpad=lp, **config["ylabel_args"])
    ax[5].set_ylabel(r"$P_\theta(t)$" +"\n"+ rf"[${Ptheta.units:P~}$]", loc=loc, labelpad=lp, **config["ylabel_args"])
    ax[6].set_ylabel(r"$P_\zeta(t)$"  +"\n"+ rf"[${Pzeta.units:P~}$]",  loc=loc, labelpad=lp, **config["ylabel_args"])
    # fmt: on

    # Zoom out Pzeta y-axis
    current_ylim = np.array(ax[6].get_ylim())
    ax[6].set_ylim(np.sort([current_ylim[0] / 3, current_ylim[1] * 3]))

    # Move yticks to the right
    for i in range(len(fig.axes)):
        ax[i].yaxis.tick_right()

    # Add common xlabel
    fig.supxlabel(f"t [{t.units:~P}]", **config["ylabel_args"])

    fig.set_tight_layout(True)
    plt.ion()
    plt.show(block=True)

    logger.info("--> Time evolutions successfully plotted.")
