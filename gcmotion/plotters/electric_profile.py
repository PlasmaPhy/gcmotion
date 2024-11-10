r"""
Electric Profile
----------------

Plots a poloidal cut of the Electric potential, as well as 2 plots of the 
electric potential and field with respect to :math:`r`.
"""

import pint
import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import electric_profile as config

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def electric_profile(
    tokamak: dict, Q: Quantity, zoom: list = [0, 1.1], units: str = "SI"
):
    r"""Plots a poloidal cut of the Electric potential, as well as 2 plots of the
    electric potential and field with respect to :math:`r`.

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
    efield = tokamak["efield"]
    Ea = efield.Ea
    EaNU = efield.EaNU
    fig = plt.figure(figsize=config["figsize"])

    # Setup x-axis lim
    psi_wallNU = (a / R) ** 2 / 2  # NU
    psi_min, psi_max = zoom[0] * psi_wallNU, zoom[1] * psi_wallNU
    aspan = abs((zoom[0] - zoom[1]) * a)

    # Grid creation and (b,i,g) calculation
    # Coordinates psi, r, theta are all in [NU]
    density = config["grid_density"]
    psi, theta = np.meshgrid(
        np.linspace(psi_min.m, psi_max.m, density),
        np.linspace(0, 2 * np.pi, density),
    )
    PhiNU = efield.PhiNU(psi, theta)
    ErNU = efield.Er(psi[0])  # 1d psi
    r = np.sqrt(2 * psi)

    # Set appropriate units:
    if units.lower() == "nu":
        R = a.to("NUmeters")
        a = a.to("NUmeters")
        aspan.ito("NUmeter")
        Phi = Q(PhiNU, "NUVolts")
        Eunits = EaNU.to_compact().units
        Er = Q(ErNU, "NUVolts/NUmeter")
    else:
        r = Q(r, "NUmeters").to("meters")
        Phi = Q(PhiNU, "NUVolts").to("Volts")
        Phiunits = np.abs(Phi[0]).max().to_compact().units
        Phi.ito(Phiunits)
        Eunits = Ea.to_compact().units
        Er = Q(ErNU, "NUVolts/NUmeter").to(Eunits)

    Ermax = Er.m.max()
    Phimax = Phi[0].m.max()

    aspectE = aspan.m / np.abs(Er).m.max()
    aspectPhi = aspan.m / np.abs(Phimax)

    # ===============================================================

    # Electric Field contour
    ax0 = fig.add_subplot(1, 3, 1, projection="polar")
    ax0.set_title(rf"Electric Potential $\Phi(r)$[{Phi.units:~P}]")
    ax0.set_rlabel_position(70)
    ax0.set_yticks([a.m], labels=[f"{a:.3g}"])
    C0 = ax0.contourf(theta, r.m, Phi.m, **config["contour_params"])
    fig.colorbar(
        C0,
        ax=ax0,
        pad=0.1,
        shrink=0.4,
        label=f"[{Phi.units:~P}]",
    )

    # Electric Field 2D plot
    ax1 = fig.add_subplot(1, 3, 2)
    ax1.set_title(rf"Electric Field $E(r)$[{Er.units:~P}]")
    ax1.set_xlabel(rf"$r$[{r.units:~P}]")
    ax1.set_ylabel(rf"$E(r)$[{Er.units:~P}]")
    ax1.plot(r[0].m, Er, **config["plot_params"])
    ax1.vlines(a, Er.min(), Ermax, **config["vline_params"])
    ax1.set_aspect(config["aspect_ratio"] * aspectE)

    # Electric Potential 2D plot
    ax2 = fig.add_subplot(1, 3, 3)
    ax2.set_title(rf"Electric Potential $\Phi(r)$[{Phi.units:~P}]")
    ax2.set_xlabel(rf"$r$[{r.units:~P}]")
    ax2.set_ylabel(rf"$\Phi(r)$[{Phi.units:~P}]")
    ax2.plot(r[0].m, Phi[0], **config["plot_params"])
    ax2.vlines(a, Phi[0].min(), Phimax, **config["vline_params"])
    ax2.set_aspect(config["aspect_ratio"] * aspectPhi)

    plt.tight_layout()
    plt.ion()
    plt.show(block=True)
