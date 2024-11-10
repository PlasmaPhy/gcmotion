r"""
Magnetic Profile
----------------

Plots a poloidal cut of the magnetic field strength and plasma current
intensities.
"""

import pint
import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import magnetic_profile as config

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def magnetic_profile(tokamak: dict, Q: Quantity, units: str = "SI"):
    r"""Plots a poloidal cut of the magnetic field strength and plasma
    current intensities.

    Parameters
    ----------
    tokamak : dict
        Dictionary containing the tokamak configuration, as passed to
        a particle.
    Q : Quantity
        The Quantity constructor.
    units : str, optional
        The units of the contour plots. Can be either "NU" or "SI".
        Defaults to "SI".
    """

    logger.info("Plotting Magnetic field profile...")
    R = tokamak["R"].to("meters")
    a = tokamak["a"].to("meters")
    B0 = tokamak["B0"].to("Tesla")
    bfield = tokamak["bfield"]

    fig, ax = plt.subplots(
        1, 3, figsize=config["figsize"], subplot_kw={"projection": "polar"}
    )

    # Grid creation and (b,i,g) calculation
    # Coordinates psi, r, theta are all in [NU]
    density = config["grid_density"]
    psi_wallNU = (a / R) ** 2 / 2  # NU
    psi, theta = np.meshgrid(
        np.linspace(0, psi_wallNU.m, density),
        np.linspace(0, 2 * np.pi, density),
    )
    b, i, g = bfield.bigNU(psi, theta)
    r = np.sqrt(2 * psi)

    # Set appropriate units:
    if units.lower() == "nu":
        R = tokamak["R"].to("NUmeters")
        a = tokamak["a"].to("NUmeters")
        B0 = tokamak["B0"].to("NUTesla")
        b = B0 * b  # NUTesla
        i = Q(i, "NUPlasma_current")
        g = Q(g, "NUPlasma_current")
    else:
        r = Q(r, "NUmeters").to("meters").m
        b = B0 * b  # Tesla
        i = Q(i, "NUPlasma_current").to("Plasma_current")
        g = Q(g, "NUPlasma_current").to("Plasma_current")

    # ===============================================================

    # Bfield contour
    ax[0].set_title(f"Magnetic field strength [{B0.units:~P}]")
    ax[0].set_rlabel_position(70)
    ax[0].set_yticks([a.m], labels=[f"{a:.3g}"])
    C0 = ax[0].contourf(theta, r, b.m, **config["bcontour_params"])
    fig.colorbar(
        C0,
        ax=ax[0],
        pad=0.1,
        shrink=0.4,
        label=f"[{B0.units:~P}]",
        ticks=np.sort([b.m.min(), B0.m, b.m.max()]),
    )

    # I contour
    ax[1].set_title(rf"Toroidal current 'I'[{i.units:~P}]")
    ax[1].set_rlabel_position(70)
    ax[1].set_yticks([a.m], labels=[f"{a:.3g}"])
    C1 = ax[1].contourf(theta, r, i.m, **config["icontour_params"])
    fig.colorbar(
        C1,
        ax=ax[1],
        pad=0.1,
        shrink=0.4,
        label=f"[{i.units:~P}]",
        ticks=np.sort([i.m.min(), i.m.mean(), i.m.max()]),
    )

    # g contour
    ax[2].set_title(rf"Poloidal current 'g'[{g.units:~P}]")
    ax[2].set_rlabel_position(70)
    ax[2].set_yticks([a.m], labels=[f"{a:.3g}"])
    C2 = ax[2].contourf(theta, r, g.m, **config["gcontour_params"])
    fig.colorbar(
        C2,
        ax=ax[2],
        pad=0.1,
        shrink=0.4,
        label=f"[{g.units:~P}]",
        ticks=np.sort([g.m.min(), g.m.mean(), g.m.max()]),
    )

    plt.tight_layout()
    plt.ion()
    plt.show(block=True)
