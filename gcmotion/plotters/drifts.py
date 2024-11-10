r"""
Plots poloidal :math:`\theta - P_\theta` and
:math:`\zeta - P_\zeta` drifts.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

Example
-------

.. code-block:: python

    gcm.drifts(cwp, lim[-np.pi, np.pi])

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger
from gcmotion.utils.plot_utils import pi_mod

from gcmotion.configuration.plot_parameters import drifts as config


def drifts(
    cwp,
    theta_lim: list = [-np.pi, np.pi],
    units: str = "SI",
    **params,
):
    r"""
    Draws 2 plots: 1] :math:`\theta-P_\theta`
    and 2] :math:`\zeta-P_\zeta`.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    theta_lim : list, optional
        Plot xlim. Must be either [0,2π] or [-π,π]. Defaults to [-π,π].
    units : str, optional
        The unit system. Can be either 'SI' or 'NU'. Defauls to "SI".
    params : dict, optional
        Extra arguements if called for many particles:

            #. plot_initial : bool, optional
                Whether or not to plot the initial point. Defaults to True
    """
    if params.get("_internal_call", False):  # dont pop it yet
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    suffix = "NU" if units == "NU" else "" if units == "SI" else ""
    logger.info(f"Plotting θ-Pθ and ζ-Pζ drifts in {"NU" if suffix=="NU" else "SI"}...")  # fmt: skip

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!
    different_colors = params.get("different_colors", False)
    plot_initial = params.get("plot_initial", True)

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    # Get all needed attributes first
    theta = getattr(cwp, "theta").copy()
    zeta = getattr(cwp, "zeta").copy()
    Ptheta = getattr(cwp, "Ptheta" + suffix).copy()
    Pzeta = getattr(cwp, "Pzeta" + suffix).copy()

    if canvas is None:
        fig, ax = plt.subplots(1, 2, **config["fig_parameters"])
        fig.tight_layout()
        canvas = (fig, ax)
        logger.debug("\tDrifts: Creating a new canvas.")
    else:  # Use external canvas
        fig, ax = canvas
        logger.debug("\tDrifts: Using existing canvas.")

    fig.suptitle(r"Drift orbits of $P_\theta - \theta$ and $P_\zeta - \zeta$")

    scatter_kw = config["scatter_args"]
    if different_colors and "color" in scatter_kw.keys():
        del scatter_kw["color"]

    # Set theta lim. Mods all thetas to 2π
    theta_min, theta_max = theta_lim
    theta_plot, ax_lim = pi_mod(theta, theta_lim)

    ax[0].scatter(theta_plot, Ptheta, **config["scatter_args"])
    ax[1].scatter(zeta, Pzeta, **config["scatter_args"])

    ax[0].set_xlabel(r"$\theta$" + f"[{theta.units:~P}]", fontsize=config["xfontsize"])  # fmt: skip
    ax[1].set_xlabel(r"$\zeta$" + f"[{zeta.units:~P}]", fontsize=config["xfontsize"])  # fmt: skip

    ax[0].set_ylabel(r"$P_\theta/\psi_{wall}$", fontsize=config["yfontsize"])  # fmt: skip
    ax[1].set_ylabel(r"$P_ζ$" + rf"$[{Ptheta.units:~#P}]$", fontsize=config["yfontsize"])  # fmt: skip

    if plot_initial:
        ax[0].scatter(theta_plot[0], Ptheta[0], c="k", s=10, zorder=3)  # fmt: skip
        ax[1].scatter(zeta[0], Pzeta[0], c="k", s=10, zorder=3)

    # Zoom out Pzeta y-axis
    if not _internal_call:
        current_ylim = np.array(ax[1].get_ylim())
        ax[1].set_ylim([current_ylim[0] / 5, current_ylim[1] * 5])

    # Set all xticks as multiples of π, and then re-set xlims (smart!)
    ticks = ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"]
    ax[0].set_xticks(np.linspace(-2 * np.pi, 2 * np.pi, 9), ticks)
    ax[0].set_xlim(ax_lim)

    # Make interactive if single particle:
    if not _internal_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    logger.info("θ-Pθ and ζ-Pζ drifts successfully plotted.")
    logger.enable("gcmotion")
