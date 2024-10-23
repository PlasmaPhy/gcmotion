r"""
Plots poloidal :math:`\theta - P_\theta` or 
:math:`\zeta - P_\zeta` drift.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

Example
-------

.. code-block:: python

    gcm.drift(cwp, angle="theta")

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils._logger_setup import logger
from gcmotion.utils.pi_mod import pi_mod

from gcmotion.configuration.plot_parameters import drift as config


def drift(
    cwp,
    angle: str = "theta",
    lim: list = [-np.pi, np.pi],
    **params,
):
    r"""Draws :math:`\theta - P_\theta` plot.

    This method is called internally by ``countour_energy()``
    as well.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    angle : str, optional
        The angle to plot. Defaults to "theta".
    lim : list, optional
        Plot xlim. Must be either [0,2π] or [-π,π]. Defaults to [-π,π].
    params : dict, optional
        Extra arguements if called for many particles:

            #. plot_initial : bool, optional
                Whether or not to plot the initial point. Defaults to True
    """
    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!
    different_colors = params.get("different_colors", False)
    plot_initial = params.get("plot_initial", True)

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    logger.info(f"Plotting {angle}-P_{angle} drift...")

    # Get all needed attributes first
    q = getattr(cwp, angle)
    P_plot = getattr(cwp, "P" + angle).copy()
    psi_wall = cwp.psi_wall

    if canvas is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        canvas = (fig, ax)
        logger.debug("\tCreating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tUsing existing canvas.")

    if angle == "theta":  # Normalize to psi_wall
        P_plot /= psi_wall

    # Set theta lim. Mods all thetas or zetas to 2π or [-π,π]
    q_plot, ax_lim = pi_mod(q, lim)

    scatter_kw = config["scatter_args"].copy()
    if different_colors:
        scatter_kw.pop("color", None)

    ax.scatter(q_plot, P_plot, **scatter_kw, zorder=2)
    ax.set_xlabel(rf"$\{angle}$", fontsize=config["xfontsize"])
    ax.set_ylabel(rf"$P_\{angle}/\psi_w$", fontsize=config["yfontsize"])

    if plot_initial:
        ax.scatter(q[0], P_plot[0], c="k", s=10, zorder=3)

    if angle == "zeta":  # Zoom out Pzeta y-axis
        current_ylim = np.array(ax.get_ylim())
        ax.set_ylim([current_ylim[0] / 5, current_ylim[1] * 5])

    # Set all xticks as multiples of π, and then re-set xlims (smart!)
    ticks = ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"]
    ax.set_xticks(np.linspace(-2 * np.pi, 2 * np.pi, 9), ticks)
    ax.set_xlim(ax_lim)

    if not _internal_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    logger.info(f"{angle}-P_{angle} drift successfully plotted.")
    logger.enable("gcmotion")
