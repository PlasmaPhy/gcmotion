import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Rectangle
from matplotlib.axes import Axes

from gcmotion.utils.logger_setup import logger
from gcmotion.configuration.plot_parameters import (
    BaseProfileContourConfig as config,
)
from gcmotion.entities.profile import Profile


def _base_profile_contour(profile: Profile, ax: Axes, **args):
    r"""Base plotting function. Only draws upon a given axis without showing
    any figures.

    Parameters
    ----------
    profile : Profile
        The profile entity.
    ax : Axes
        The ax upon which to draw.
    args : dict
        The optional arguement dictionary.
    """

    logger.info("==> Plotting Base Profile Contour...")

    # Unpack parameters
    # Default values are also defined here
    _thetalim = args.get("thetalim")
    _psilim = args.get("psilim")
    levels = args.get("levels")
    flux_units = args.get("flux_units")
    E_units = args.get("E_units")
    potential = args.get("potential")
    wall = args.get("wall")
    cursor = args.get("cursor")

    # Setup meshgrid
    # The Energy calculation and y axis limits are always set with respect to
    # psilim. Pthetalim only acts as a way to move the axis window at the end,
    # or zoom to a certain area with respect to Ptheta.
    thetalim = profile.Q(_thetalim, "radians")
    psilim = profile.Q(_psilim, "psi_wall").to(flux_units)
    psigrid, thetagrid = np.meshgrid(
        np.linspace(psilim[0], psilim[1], config.grid_density),
        np.linspace(thetalim[0], thetalim[1], config.grid_density),
    )

    # Calculate Energy values
    Energy = profile.findEnergy(
        psigrid, thetagrid.m, E_units, potential=potential
    )

    # Define data to be plotted
    # Take only the magnitudes to supress warnings
    data = {
        "theta": thetagrid.m,
        "flux": psigrid.m,
        "Energy": Energy.m,
    }

    # =============================
    # Contour related configuration
    # =============================

    # Locator setup
    locator = (
        ticker.LogLocator(base=config.log_base, numticks=levels)
        if config.locator == "log"
        else ticker.MaxNLocator(nbins=levels)
    )

    kw = {
        "cmap": config.cmap,
        "locator": locator,
        "zorder": 0,
    }

    # Contour plot
    C = ax.contourf("theta", "flux", "Energy", data=data, **kw)
    plt.rcParams["axes.autolimit_mode"] = "round_numbers"
    # Setup labels.
    # Also add a second axis for Ptheta
    ax.set_ylabel(rf"$\psi$ [{psigrid.units:.4g~P}]", size=config.labelsize)
    ax.set_xlabel(r"$\theta$ [radians]", size=config.labelsize)
    ax.tick_params(labelsize=config.ticksize)
    ax.set_xticks(
        np.linspace(-2 * np.pi, 2 * np.pi, 9),
        ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"],
        size=config.ticksize,
    )
    ax.set_xlim(thetalim.m)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(config.ticknum))

    ax2 = ax.twinx()
    ax2.set_ylabel(
        rf"$P_\theta$ [{psigrid.units:.4g~P}]", size=config.labelsize
    )
    psiticks = profile.Q(ax.get_yticks(), flux_units)
    Pthetamin = profile.findPtheta(psiticks.min())
    Pthetamax = profile.findPtheta(psiticks.max())
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(config.ticknum))
    ax2.set_ylim([Pthetamin.m, Pthetamax.m])
    ax2.tick_params(labelsize=config.ticksize)

    # Add a shade above psi_wall
    if wall:
        rect = Rectangle(
            (thetalim.min().m, profile.Q(1, "psi_wall").to(flux_units).m),
            width=abs(thetalim[1].m - thetalim[0].m),
            height=psigrid[0][-1].m,
            alpha=0.2,
            color="k",
        )
        ax.add_patch(rect)

    # Re-format cursor
    # if cursor:
    #     Pthetagrid = profile.findPtheta(psigrid)
    #
    #     ax.format_coord = grid_cursor_format_coord(
    #         thetagrid.m,
    #         psigrid.m,
    #         Energy.m,
    #         profile,
    #         flux_units,
    #         E_units,
    #     )

    # Return the contour object
    return C


# def grid_cursor_format_coord(
#     theta: np.ndarray,
#     psi: np.ndarray,
#     Energy: np.ndarray,
#     profile: Profile,
#     flux_units: str,
#     z_label: str = "",
# ):
#     r"""Creates a data cursor on the contour plot."""
#
#     thetaFlat, psiFlat, EnergyFlat = (
#         theta.flatten(),
#         psi.flatten(),
#         Energy.flatten(),
#     )
#
#     def fmt(x, y):
#         # get closest point with known data
#         dist = np.linalg.norm(np.vstack([thetaFlat - x, psiFlat - y]), axis=0)
#         idx = np.argmin(dist)
#         # y2 = PthetaFlat[idx]
#         y2 = profile.findPtheta(profile.Q(y, flux_units))
#         z = EnergyFlat[idx]
#         return (
#             f"theta={x:.3g},"
#             + f"psi={y:.4g},Ptheta={y2:.4g},"
#             + f"E={z:.4g}[{z_label}]"
#         )
#
#     return fmt
