r"""
Plots poloidal :math:`\theta - P_\theta` and
:math:`\zeta - P_\zeta` drifts of a collection of particles.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

Example
-------

.. code-block:: python

    params = {"different_colors": False, "plot_initial": True}
    gcm.collection_drift(
        collection, angle = "theta",
        lim = [-np.pi, np.pi], params=params
    )

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.plotters.drifts import drifts

from gcmotion.utils.logger_setup import logger


def collection_drifts(
    collection,
    theta_lim: list = [-np.pi, np.pi],
    **params,
):
    r"""
    Plots poloidal :math:`\theta - P_\theta` and
    :math:`\zeta - P_\zeta` drifts of a collection of particles.

    Parameters
    ----------
    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles
    theta_lim : list, optional
        Plot xlim. Must be either [0,2π] or [-π,π]. Defaults to [-π,π].
    params : dict, optional
        Extra plotting parameters:

            * different_colors : bool
                Whether or not not use different colors for every drift.
                Defaults to False.
            * plot_initial: bool
                Whether or not to plot the starting points of each drift.
                Defaults to True.
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
        fig, ax = plt.subplots(1, 2, figsize=(12, 6))
        fig.tight_layout()
        canvas = (fig, ax)
        logger.debug("\tCreating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tUsing existing canvas.")

    for p in collection.particles:
        drifts(
            p,
            theta_lim=theta_lim,
            _internal_call=True,
            canvas=canvas,
            **params,
        )

    if not _internal_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    logger.enable("gcmotion")
    logger.info("--> Collection drift successfully plotted.")
