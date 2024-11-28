import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Rectangle

from gcmotion.utils.logger_setup import logger
from gcmotion.configuration.plot_parameters import (
    BaseProfileContourConfig as config,
)
from gcmotion.entities.profile import Profile


def base_profile_contour(profile: Profile, **args):
    r"""Plots a basic theta-psi/Ptheta-Hamiltonian contour.

    Parameters
    ----------
    thetalim : list, optional
        The :math:`\theta` span in radians. Defaults to [-π,π].
    psilim : list, optional
        The :math:`\psi` span in units of *psi_wall*. Defaults to [0, 1.2].
    levels : int, optional
        The number of contour lines. Defaults to 25.
    flux_units : str, optional
        The units of the psi/Ptheta axis. Can be "Magnetic_flux"(same as "Tesla
        * meters^2"), "NUmagnetic_flux" (same as "NUTesla * NUmeters^2"),
        "psi_wall", "NUpsi_wall", or any other pint unit with the same
        dimensionality. Defaults to "Tesla * meters^2".
    E_units : str, optional
        The Energy units. Can be "eV", "keV", "Joule", "NUJoule" or any other
        pint unit with the same dimensionality. Defaults to "keV".
    potential : bool, optional
        Whether or not to add the potential Φ term when calculating the Energy.
        Defaults to True.
    wall : bool, optional
        Whether or not to shade the area above :math:`\psi_wall`. Defaults to
        True.

    Notes
    -----

    This function cal also take some extra, hidden parameters to add
    functionality when used inside a plotting pipeline:

        #.  canvas
        #.  show
    """

    logger.info("==> Plotting Base Profile Contour...")

    # Unpack parameters
    # Default values are also defined here
    _thetalim = args.get("thetalim", config.thetalim)
    _psilim = args.get("psilim", config.psilim)
    levels = args.get("levels", config.levels)
    flux_units = args.get("flux_units", config.flux_units)
    E_units = args.get("units", config.E_units)
    potential = args.get("potential", config.potential)
    wall = args.get("wall", config.wall)
    show = args.get("show", True)
    canvas = args.get("canvas", None)
    del args

    # Create figure
    # An ax is created unless an already existing one is passed in args.
    figkw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
    }
    fig = canvas[0] if canvas is not None else plt.figure(**figkw)
    ax = canvas[1] if canvas is not None else fig.subplots()

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
        psigrid, thetagrid.m, E_units, potential=potential)

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
        "zorder":  0,
    }

    # Contour plot
    ax.contourf("theta", "flux", "Energy", data=data, **kw)
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    # Setup labels.
    # Also add a second axis for Ptheta
    ax.set_ylabel(rf"$\psi$ [{psigrid.units:.4g~P}]", size=config.labelsize)
    ax.set_xlabel(r"$\theta$ [radians]", size=config.labelsize)
    ax.tick_params(labelsize=config.ticksize)
    ax.set_xticks(
        np.linspace(-2*np.pi, 2*np.pi, 9),
        ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"],
        size=config.ticksize
    )
    ax.set_xlim(thetalim.m)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(config.ticknum))

    ax2 = ax.twinx()
    ax2.set_ylabel(rf"$P_\theta$ [{psigrid.units:.4g~P}]",
                   size=config.labelsize)
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

    if show:
        plt.show()
