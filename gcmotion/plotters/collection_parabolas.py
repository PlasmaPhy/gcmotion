r"""
In the absence of electric field and in LAR configuration, it plots
the orbit type parabolas as described by White.

Also plots the particles' orbit points.

Example
-------

.. code-block:: python

    gcm.collection_parabolas(collection)

.. rubric:: Function:
    :heading-level: 4

"""

import matplotlib.pyplot as plt

from gcmotion.plotters.parabolas import parabolas, orbit_point

from gcmotion.utils.logger_setup import logger


def collection_parabolas(
    collection,
    **params,
):
    r"""
    Plots the poloidal :math:`\theta - P_\theta` drifts of
    a collection of particles.

    Parameters
    ----------
    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles
    params : dict, optional
        Extra plotting parameters:

            * different_colors : bool
                Whether or not not use different colors for every drift.
                Defaults to False.
            * labels : bool
                Whether or not to plot a label above particle.
            * autoscale : bool
                Whether to autoscale the x axis in order to get
                all the points int view, or limit the plot to the
                parabolas' area.
    """
    if collection.has_efield and collection.is_lar:
        string = (
            "At least 1 particle inside the collection has efield. Aborting..."
        )
        print(string)
        logger.warning(string)
        return

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!
    plot_label = params.get("plot_label", True)
    autoscale = params.get("autoscale", False)

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")
    logger.info("Plotting Collection drift...")

    if canvas is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        canvas = (fig, ax)
        logger.debug("\tCreating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tUsing existing canvas.")

    for p in collection.particles:
        orbit_point(
            cwp=p,
            _internal_call=True,
            canvas=canvas,
            plot_label=plot_label,
        )

    parabolas(
        collection.particles[0],
        plot_point=False,
        _internal_call=True,
        canvas=canvas,
        autoscale=autoscale,
    )

    if not _internal_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    logger.enable("gcmotion")
    logger.info("--> Collection drift successfully plotted.")
