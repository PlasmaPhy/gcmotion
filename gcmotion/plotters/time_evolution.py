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

    gcm.time_evolution(cwp, percentage=20, units="s")

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils._logger_setup import logger

from gcmotion.configuration.plot_parameters import time_evolution as config


def time_evolution(cwp, percentage: int = 100, units: str = "s"):
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
        The time units. Can be either 's' for seconds or 'NU' for
        normalized units. Defauls to "s".
    """
    logger.info("Plotting time evolutions...")

    # Get all needed attributes first (.copy() if needed!)
    t = cwp.t_eval.copy()
    theta = cwp.theta
    psi = cwp.psi
    zeta = cwp.zeta
    rho = cwp.rho
    psip = cwp.psip
    Ptheta = cwp.Ptheta
    Pzeta = cwp.Pzeta
    psip_wall = cwp.psip_wall

    if percentage < 1 or percentage > 100:
        percentage = 100
        print("Invalid percentage. Plotting the whole thing.")
        logger.warning("Invalid percentage: Plotting the whole thing...")

    points = int(np.floor(t.shape[0] * percentage / 100) - 1)

    if units == "s":
        t /= cwp.w0

    # Plotting
    fig, ax = plt.subplots(7, 1, **config["fig_parameters"])
    ax[0].set_title("Time evolution of dynamical variables", c="b")
    ax[5].set_title("Time evolution of canonical momenta", c="b")

    ax[0].scatter(t[:points], theta[:points], **config["scatter_args"])
    ax[1].scatter(t[:points], zeta[:points], **config["scatter_args"])
    ax[2].scatter(t[:points], psi[:points], **config["scatter_args"])
    ax[3].scatter(t[:points], psip[:points], **config["scatter_args"])
    ax[4].scatter(t[:points], rho[:points], **config["scatter_args"])
    ax[5].scatter(t[:points], Ptheta[:points], **config["scatter_args"])
    ax[6].scatter(t[:points], Pzeta[:points], **config["scatter_args"])

    ax[0].set_ylabel(r"$\theta(t)$", **config["ylabel_args"])
    ax[1].set_ylabel(r"$\zeta(t)$", **config["ylabel_args"])
    ax[2].set_ylabel(r"$\psi(t)$", **config["ylabel_args"])
    ax[3].set_ylabel(r"$\psi_p(t)$", **config["ylabel_args"])
    ax[4].set_ylabel(r"$\rho(t)$", **config["ylabel_args"])
    ax[5].set_ylabel(r"$P_\theta(t)\quad$", **config["ylabel_args"])
    ax[6].set_ylabel(r"$P_\zeta(t)$", **config["ylabel_args"])
    ax[6].set_ylim([-psip_wall, psip_wall])

    if units == "NU":
        fig.supxlabel("t [NU]", **config["ylabel_args"])
    elif units == "s":
        fig.supxlabel("t [s]", **config["ylabel_args"])

    fig.set_tight_layout(True)
    plt.ion()
    plt.show(block=True)

    logger.info("--> Time evolutions successfully plotted.")
