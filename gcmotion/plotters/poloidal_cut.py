r"""
Plots the Poloidal view of the torus and the 
:math:`\theta - P_\theta` drift.

Example
-------

.. code-block:: python

    gcm.poloidal_cut(cwp, plot_initial=False)
    
.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap, colors

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import poloidal_cut as config

from gcmotion.utils.canonical_to_toroidal import canonical_to_toroidal


def poloidal_cut(
    cwp, wall_shade: bool = True, plot_axis: bool = True, **params
):
    r"""Plots the poloidal cut of the orbit in polar coordinates.

    :meta public:

    Parameters
    -----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    wall_shade : bool, optional
        Whether or not to shade the area  close to:math:`r>r_{wall}`.
        Defaults to True
    plot_axis : bool, optional
        Whether or not to plot the magnetic axis. Defaults to True.
    params : dict, optional
        Extra plotting parameters:

            * different_colors : bool
                Whether or not not use different colors for every drift.
                Defaults to False.
            * plot_initial : bool
                Whether or not to plot the starting points of each drift.
                Defaults to True.
            * plot_axis : bool
                Whether or not to plot the magnetic axis
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

    if canvas is None:
        fig, ax = plt.subplots(
            figsize=(12, 12),
            subplot_kw={"projection": "polar"},
        )
        canvas = (fig, ax)
        logger.debug("\tCreating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tUsing existing canvas.")

    logger.info("Plotting 2D torus sections...")
    # Configure torus dimensions and orbit and store internally
    Rtorus, atorus, r_torus, theta_torus, z_torus = canonical_to_toroidal(
        cwp, percentage=100, truescale=True
    )

    r_plot1 = r_torus

    orbit_kw = config["orbit_kw"].copy()
    if different_colors:
        orbit_kw.pop("color", None)

    if plot_initial:
        ax.scatter(theta_torus[0], r_plot1[0], c="k", s=10, zorder=3)

    if plot_axis:
        s = config["axis_size"]
        ax.scatter(0, 0, s=s, **config["axis_kwargs"])

    # Orbits
    ax.scatter(theta_torus, r_plot1, **orbit_kw, zorder=-1)

    ax.set_ylim(bottom=0, top=atorus * 1.1)
    ax.grid(False)
    ax.set_title("Poloidal View", c="b")
    ax.set_xlabel(r"$\sqrt{2\psi} - \theta$")
    ax.tick_params(labelsize=8)
    ax.set_axis_off()

    logger.info("--> 2D torus sections plotted successfully.")

    # Hard y limit
    top = plt.gca().get_ylim()[1]
    plt.autoscale(axis="y")
    if top > atorus * 2:
        plt.ylim(top=atorus * 2)

    if not _internal_call:
        _wall(canvas, atorus, wall_shade)  # Wall
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)


def _wall(canvas, atorus, wall_shade):
    """Creates the wall and shade of the torus.

    :meta private:

    Parameters
    ----------
    canvas : 2-tuple
        2-tuple containing (fig, ax)
    atorus : float
        The torus' minor radius in [m]
    """
    fig, ax = canvas

    wall_points = config["wall_points"]

    # Torus Walls
    ax.scatter(
        np.linspace(0, 2 * np.pi, wall_points),
        atorus * np.ones(wall_points),
        **config["wall_kw"],
    )

    # Shade
    if wall_shade:
        x = np.linspace(0, 2 * np.pi, 100)
        y = np.linspace(0.8 * atorus, atorus, 100)

        cmap = get_cmap("Greys")
        norm = colors.Normalize(vmin=0.8 * atorus, vmax=atorus)

        for i in range(len(x) - 1):
            plt.fill_between(
                x,
                y[i],
                y[i + 1],
                color=cmap(norm(y[i])),
                alpha=0.1,
                zorder=-1,
            )
