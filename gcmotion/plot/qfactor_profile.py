r"""
Afactor Profile
----------------

Plots :math:`q(\psi)` and :math:`\psi_p(\psi)`.
"""

import pint
import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import qfactor_profile as config

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def qfactor_profile(
    tokamak: dict, Q: Quantity, zoom: list = [0, 1.1], units: str = "SI"
):
    r"""Plots :math:`q(\psi)` and :math:`\psi_p(\psi)`.

    Parameters
    ----------
        tokamak : dict
            Dictionary containing the tokamak configuration, as passed to
            a particle.
        Q : Quantity
            The Quantity constructor.
        zoom: list
            The x-axis zoom limits.
        units : str, optional
            The units of the contour plots. Can be either "NU" or "SI".
            Defaults to "SI".
    """

    logger.info("Plotting Magnetic field profile...")
    R = tokamak["R"].to("meters")
    a = tokamak["a"].to("meters")
    qfactor = tokamak["qfactor"]

    fig = plt.figure(figsize=config["figsize"])

    # Setup x-axis lim
    psi_wall = (a.m / R.m) ** 2 / 2
    psi_wall = Q(psi_wall, "NUmagnetic_flux")  # NU
    psi_min, psi_max = zoom[0] * psi_wall, zoom[1] * psi_wall

    # Calculate values and units
    psi = np.linspace(psi_min.m, psi_max.m, 500)
    q = qfactor.solverqNU(psi)
    psip = qfactor.psipNU(psi)

    psi = Q(psi, "NUMagnetic_flux")
    psip = Q(psip, "NUMagnetic_flux")

    # Set appropriate units:
    if units.lower() == "nu":
        R.ito("NUmeters")
        a.ito("NUmeters")
    else:
        psi_wall.ito("Magnetic_flux")
        psi.ito("Magnetic_flux")
        psip.ito("Magnetic_flux")

    # ===============================================================

    # q(ψ) 2D plot
    ax0 = fig.add_subplot(1, 2, 1)
    ax0.set_title(r"$q(\psi)$")
    ax0.set_xlabel(rf"$\psi$[{psi.units:~P}]")
    ax0.set_ylabel(r"$q(\psi)$")
    ax0.plot(psi.m, q, **config["plot_params"])
    ax0.vlines(psi_wall, q.min(), q.max(), **config["vline_params"])

    # Electric Potential 2D plot
    ax1 = fig.add_subplot(1, 2, 2)
    ax1.set_title(r"$\psi_p(\psi)$")
    ax1.set_xlabel(rf"$\psi$[{psi.units:~P}]")
    ax1.set_ylabel(rf"$\psi_p$[{psip.units:~P}]")
    ax1.plot(psi.m, psip.m, **config["plot_params"])
    ax1.vlines(psi_wall, psip.min(), psip.max(), **config["vline_params"])

    plt.tight_layout()
    plt.ion()
    plt.show(block=True)
