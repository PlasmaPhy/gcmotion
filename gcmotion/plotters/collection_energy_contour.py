r"""
Plots the energy contour lines of the :math:`\theta-P_\theta`
plain of the Hamiltonian of a collection of particles with the same
required constants of motion.

Can also plot the current particle's :math:`\theta-P_\theta` drift.

We can optionally remove the :math:`\Phi` term from the Hamiltonian,
to observe the change of the drifts with the added electric field.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

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
from gcmotion.plotters.energy_contour import _cbar

from gcmotion.configuration.plot_parameters import energy_contour as config

from gcmotion.utils.logger_setup import logger


def collection_energy_contour(
    collection,
    theta_lim: list = [-np.pi, np.pi],
    Ptheta_lim: str | list = "auto",
    plot_drift: bool = True,
    contour_Phi: bool = True,
    units: str = "keV",
    levels: int = None,
    wall_shade: bool = True,
    **params,
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
    Ptheta_lim : list | str, optional
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
        The energy units. Must be 'NU', 'eV' or 'keV'. Defaults to "keV".
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
            * plot_initial : bool
                Whether or not to plot the starting points of each drift.
                Defaults to True.
            * adjust_cbar : bool
                Whether or not to abjust color bar ticks and boundaries to
                the particles' energies.
    """
    suffix = "NU" if units == "NU" else "" if units == "SI" else ""

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.get("canvas", None)  # POP!
    different_colors = params.get("different_colors", False)
    plot_initial = params.get("plot_initial", True)
    adjust_cbar = params.get("adjust_cbar", False)

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    def params_ok() -> bool:
        """Checks for the validity of the parameters.

        Returns
        -------
        bool
            The check result
        """
        must_be_the_same = [
            "R",
            "a",
            "qfactor",
            "bfield",
            "efield",
            "species",
            "mu/muB",
            "Pzeta0",
        ]
        can_be_different = [
            key
            for key in collection.params.keys()
            if key not in must_be_the_same
        ]

        for quantity in must_be_the_same:
            if getattr(collection, "multiple_" + quantity):
                error_str = f"Only the variables {can_be_different} may vary from particle to particle."
                print(error_str)
                logger.error(error_str)
                return False

        return True

    def plot():
        """Does the actual plotting"""

        logger.info("Plotting Collection energy contour...")

        nonlocal Ptheta_lim, _internal_call, canvas, different_colors, plot_initial

        if canvas is None:
            fig = plt.figure(**config["fig_parameters"])
            ax = fig.add_subplot(111)
            canvas = (fig, ax)
            logger.debug("\tCreating a new canvas.")
        else:
            fig, ax = canvas
            logger.debug("\tUsing existing canvas.")

        if plot_drift:
            collection_drift(
                collection,
                angle="theta",
                theta_lim=theta_lim,
                units=units,
                _internal_call=True,
                canvas=canvas,
                **params,
            )

        # if Ptheta_lim == "auto":
        #     minpsiNU = min(
        #         getattr(collection[_], "psiNU").min()
        #         for _ in range(collection.n)
        #     )
        #     maxpsiNU = max(
        #         getattr(collection[_], "psiNU").max()
        #         for _ in range(collection.n)
        #     )
        #     psiNUlim = [minpsiNU, maxpsiNU]
        # print(psiNUlim)
        logger.info("\t\tPlotting single energy contour...")
        logger.disable("gcmotion")
        C = energy_contour(
            collection.particles[0],
            theta_lim=theta_lim,
            Ptheta_lim=Ptheta_lim,
            plot_drift=False,
            contour_Phi=contour_Phi,
            units=units,
            levels=levels,
            wall_shade=wall_shade,
            _internal_call=True,
            canvas=canvas,
            **params,
        )
        logger.enable("gcmotion")
        logger.info("\t-->Single energy contour successfully plotted.")

        # Plot colorbar
        Esuffix = "NU" if suffix == "NU" else "keV"
        energies = np.array(
            [
                getattr(collection.particles[_], "E" + Esuffix).magnitude
                for _ in range(collection.n)
            ]
        )
        cbar = _cbar(energies=energies, units=units, canvas=canvas, C=C)

        if adjust_cbar:
            boundaries = [min(energies), max(energies)]
            cbar.ax.set_yticks(energies)
            cbar.ax.set_ylim(boundaries)
            cbar.extendfrac = "auto"

        logger.debug("\tCollection call. Adding energy labels...")

        logger.info("--> Collection energy contour plotted successfully.")

        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)

    if params_ok():
        plot()
