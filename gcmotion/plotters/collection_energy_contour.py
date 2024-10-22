r"""
Plots the energy contour lines of the :math:`\theta-P_\theta`
plain of the Hamiltonian of a collection of particles with the same
required constants of motion.

Can also plot the current particle's :math:`\theta-P_\theta` drift.

We can optionally remove the :math:`\Phi` term from the Hamiltonian,
to observe the change of the drifts with the added electric field.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

The optional arguements are only used when plotting drifts
from multiple particles in the same canvas.

Example
-------

.. code-block:: python

    gcm.energy_contour(
        cwp, theta_lim = [-np.pi ,np.pi], psi_lim="auto", 
        plot_drift=True, contour_Phi=True, units="keV", 
        levels=20
    )

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.plotters.collection_drift import collection_drift
from gcmotion.plotters.energy_contour import energy_contour
from gcmotion.plotters.energy_contour import _cbar_energy

from gcmotion.configuration.plot_parameters import energy_contour as config

from gcmotion.utils._logger_setup import logger, add_logger


def collection_energy_contour(
    collection,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: str | list = "auto",
    plot_drift: bool = True,
    contour_Phi: bool = True,
    units: str = "keV",
    levels: int = None,
    wall_shade: bool = True,
    params: dict = {},
):
    r"""
    Plots the energy contour lines of the :math:`\theta-P_\theta`
    plain of the Hamiltonian of a collection of particles with the same
    required constants of motion.

    Can also plot the current particle's :math:`\theta-P_\theta` drift.

    :meta public:

    Parameters
    ----------

    theta_lim : list, optional
        Plot xlim. Must be either [0,2π] or [-π,π]. Defaults to [-π,π].
    psi_lim : list | str, optional
        If a list is passed, it plots between the 2 values relative to
        :math:`\psi_{wall}`. If "auto" is passed, it automatically sets
        the optimal :math:`\psi` limits. Defaults to 'auto'.
    plot_drift : bool, optional
        Whether or not to plot :math:`\theta-P_\theta` drift on top.
        Defaults to True.
    contour_Phi : bool, optional
        Whether or not to add the Φ term in the energy contour.
        Defaults to True.
    units : str, optional
        The energy units. Must be 'normal', 'eV' or 'keV'. Defaults to "keV".
    levels : int, optional
        The number of contour levels. Defaults to Config setting.
    wall_shade : bool, optional
        Whether to shade the region :math:`\psi/\psi_{wall} > 1`.
        Defaults to True.
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
    """

    def params_ok() -> bool:
        """checks for the validity of the parameters

        Returns
        -------
        bool
            The check result
        """
        must_be_the_same = ["R", "a", "q", "Bfield", "Efield", "species", "mu", "Pz0"]
        can_be_different = [key for key in collection.params.keys() if key not in must_be_the_same]

        if not collection._check_multiples(must_be_the_same):
            error_str = f"Only the variables {can_be_different} may vary from particle to particle."
            print(error_str)
            logger.error(error_str)
            return False
        return True

    def plot():
        """Does the actual plotting"""

        logger.info("Plotting Collection energy contour...")

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        canvas = (fig, ax)
        params["canvas"] = canvas

        # plt.figure()

        logger.info("\t\tPlotting single energy contour...")
        logger.remove()
        C = energy_contour(
            collection.particles[int(len(collection.particles) / 2)],
            theta_lim=theta_lim,
            psi_lim=psi_lim,
            plot_drift=plot_drift,
            contour_Phi=contour_Phi,
            units=units,
            levels=levels,
            wall_shade=wall_shade,
            params=params,
        )
        add_logger()
        logger.info("\t-->Single energy contour successfully plotted.")

        if plot_drift:
            collection_drift(collection, "theta", theta_lim, params=params, external_call=True)

        # Colorbar and labels
        cbar = fig.colorbar(C, ax=params["canvas"][1], fraction=0.03, pad=0.1, label=f"E[{units}]")
        cbar_kw = {
            "linestyle": "-",
            "zorder": 3,
            "color": config["cbar_color"],
        }
        for p in collection.particles:
            E_cbar = _cbar_energy(p, units)
            cbar.ax.plot([0, 1], [E_cbar, E_cbar], **cbar_kw)
        logger.debug("\tCollection call. Adding energy labels...")

        logger.info("--> Collection energy contour plotted successfully.")

        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    if params_ok():
        plot()
