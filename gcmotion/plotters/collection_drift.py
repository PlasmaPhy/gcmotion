r"""
Plots the poloidal :math:`\theta - P_\theta`  or 
:math:`\zeta - P_\zeta` drifts of a collection of particles.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

The optional arguements are only used when called internally
from other functions.

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

from gcmotion.plotters.drift import drift

from gcmotion.utils._logger_setup import logger


def collection_drift(
    collection,
    angle: str = "theta",
    lim: list = [-np.pi, np.pi],
    params={},
    external_call: bool = False,
):
    r"""
    Plots the poloidal :math:`\theta - P_\theta` drifts of
    a collection of particles.

    Parameters
    ----------

    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles
    angle : str, otpional
        The angle to be plotted. Defaults to "theta".
    lim : list, optional
        The x-axis (angle) limits. Defaults to [-π,π].
    params : dict, optional
        Extra plotting parameters:

            * different_colors : bool
                Whether or not not use different colors for every drift.
                Defaults to False.
            * plot_initial: bool
                Whether or not to plot the starting points of each drift.
                Defaults to True.
            * canvas: 2-tuple:
                Whether or not to plot upon an existing (fig, ax) canvas.
                Usually only used internally. Defaults to None, which
                creates a new canvas.
    external_call : bool
        Flag for the plotter to know if its called on its own or from another
        plotter. Leave as is, only meant to be used internally.
    """
    logger.info("Plotting Collection drift...")

    # Extr Parameters unpacking
    params["canvas"] = params.get("canvas", None)
    params["different_colors"] = params.get("different_colors", False)
    params["plot_initial"] = params.get("plot_initial", True)

    if params["canvas"] is None:
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        params["canvas"] = (fig, ax)
    else:  # Use external canvas
        fig, ax = params["canvas"]

    for p in collection.particles:
        logger.disable("gcmotion")
        drift(p, angle, lim, params=params)

    logger.enable()()
    logger.info("--> Collection drift successfully plotted.")

    if not external_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)
