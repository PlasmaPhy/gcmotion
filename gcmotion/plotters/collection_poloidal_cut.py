r"""
Plots the Poloidal view of the torus and the 
:math:`\theta - P_\theta` drift.

Example
-------

.. code-block:: python

    gcm.collection_poloidal_cut(cwp)
    
.. rubric:: Function:
    :heading-level: 4

"""

import matplotlib.pyplot as plt

from gcmotion.plotters.poloidal_cut import poloidal_cut, _wall

from gcmotion.utils.logger_setup import logger


def collection_poloidal_cut(
    collection,
    wall_shade: bool = True,
    plot_axis: bool = True,
    **params,
):
    r"""
    Plots poloidal :math:`\theta - P_\theta` and
    :math:`\zeta - P_\zeta` drifts of a collection of particles.

    Parameters
    ----------
    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles
    wall_shade : bool, optional
        Whether or not to shade the area :math:`r>r_\wall`.
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

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")
    logger.info("Plotting Collection drift...")

    if canvas is None:
        fig, ax = plt.subplots(
            figsize=(12, 12), subplot_kw={"projection": "polar"}
        )
        fig.tight_layout()
        canvas = (fig, ax)
        logger.debug("\tCreating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tUsing existing canvas.")

    atorus = collection.particles[0].a  # should be the same, obviously
    _wall(canvas=canvas, atorus=atorus, wall_shade=wall_shade)

    for p in collection.particles:
        poloidal_cut(
            cwp=p,
            wall_shade=wall_shade,  # ignored
            plot_axis=plot_axis,
            _internal_call=True,
            canvas=canvas,
            **params,
        )

    # Hard y limit
    if not _internal_call:
        top = plt.gca().get_ylim()[1]
        plt.autoscale(axis="y")
        if top > atorus * 2:
            plt.ylim(top=atorus * 2)

    if not _internal_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    logger.enable("gcmotion")
    logger.info("--> Collection drift successfully plotted.")
