import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import ticker  # for contour plot locator

from gcmotion.utils.logger_setup import logger
from gcmotion.utils.energy_Ptheta import energy_Ptheta

from ._tools import grid_cursor_format_coord

from gcmotion.configuration.plot_parameters import ProfileContourConfig


def profile_contour(
    profile,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: list = [0, 1.2],
    contour_Phi: bool = True,
    units: str = "SI",
    levels: int = None,
    wall_shade: bool = True,
    _final: bool = True,
):
    r"""
    Plots the energy contour lines of the :math:`\theta-P_\theta`
    plain of the Hamiltonian.

    :meta public:

    Parameters
    ----------
    profile: :py:class:`~gcmotion.classes.profile.Profile`
        The Profile object.
    theta_lim : list, optional
        x-axis span. Must be either [0,2π] or [-π,π]. Defaults to [-π,π].
    psi_lim : list, optional
        List containing the :math:`\psi` limits with respect to
        :math:`\psi_{wall}`. Defaults to [0,1.2]
    contour_Phi : bool, optional
        Whether or not to add the Φ term in the energy contour.
        Defaults to True.
    units : str, optional
        The unit system. Can be either 'SI' or 'NU'. Defauls to "SI".
    levels : int, optional
        The number of contour levels. Defaults to Config setting.
    wall_shade : bool, optional
        Whether to shade the region :math:`\psi/\psi_{wall} > 1`.
        Note that this only makes sense in the case of :math:`I=0`, where
        :math:`P_\theta=\psi`. Defaults to True.
    _final : bool, optional
        Whether or not this plotter is the final one in the pipeline, and
        plt.show() can execute. Pass it as "False" if you want to chain many
        plotters to draw in the same canvas.
    """

    suffix = "NU" if units == "NU" else ""
    logger.info(
        f"==> Plotting profile contour plot in {
                "NU" if suffix == "NU" else "SI"}..."
    )

    # Generate configuration and overwrite attributes if needed
    config = ProfileContourConfig()
    config.levels = levels if isinstance(levels, int) else config.levels

    if config.locator == "log":
        locator = ticker.LogLocator(
            base=config.log_base, numticks=config.levels
        )
    else:
        locator = ticker.MaxNLocator(nbins=config.levels)

    # Get all needed profile attributes
    Q = getattr(profile, "Q")
    psi_wallNU = getattr(profile, "psi_wallNU").m

    # Create Canvas
    fig_kw = dict(config.fig_kw)
    fig = plt.figure(**fig_kw)
    ax = fig.add_subplot(111)
    canvas = (fig, ax)
    logger.debug("\tProfile Contour: Creating a new canvas.")

    # Set psi_lim to NU units with respect to psi_wall
    psi_limNU = psi_lim[0] * psi_wallNU, psi_lim[1] * psi_wallNU

    # Calculate Energy AND Ptheta grid values in [NU] based on the
    # psi limits calculated above.
    grid_density = config.grid_density
    theta_grid, psiNU_grid = np.meshgrid(
        np.linspace(theta_lim[0], theta_lim[1], grid_density),
        np.linspace(psi_limNU[0], psi_limNU[1], grid_density),
    )

    WNU, PthetaNU = energy_Ptheta(
        theta=theta_grid,
        psi=psiNU_grid,
        profile=profile,
        contour_Phi=contour_Phi,
    )

    # Quantify them, convert them to SI if needed and calculate their spans
    WNU = Q(WNU, "NUJoule")
    PthetaNU = Q(PthetaNU, "NUMagnetic_flux")

    if units == "NU":
        E_units = "NU"
        W_grid = WNU.m
        Ptheta_grid = PthetaNU.m
    elif units == "SI":
        E_units = "keV"
        W_grid = WNU.to("kev").m
        Ptheta_grid = PthetaNU.to("Magnetic_flux").m

    Wspan = np.array([W_grid.min(), W_grid.max()])
    Pthetaspan = [Ptheta_grid.min(), Ptheta_grid.max()]

    # Contour plot
    contour_kw = {
        "vmin": Wspan[0],
        "vmax": Wspan[1],
        "levels": config.levels,
        "cmap": config.cmap,
        "locator": locator,
        "zorder": 1,
    }

    # Make sure Ptheta is also in the required units
    C = ax.contourf(theta_grid, Ptheta_grid, W_grid, **contour_kw)
    ax.set_xlabel(r"$\theta$ [rads]")
    ax.set_ylabel(
        rf"$P_\theta$ [{Q(suffix+"Magnetic_flux").units:~P}]", rotation=90
    )
    ticks = ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"]
    plt.xticks(np.linspace(-2 * np.pi, 2 * np.pi, 9), ticks)
    ax.set(xlim=theta_lim, ylim=Pthetaspan)
    ax.set_facecolor("white")

    # Setup colorbar
    Cbar = fig.colorbar(
        C,
        ax=ax,
        fraction=0.03,
        pad=0.1,
        label=f"E[{E_units}]",
        format="{x:.4g}",
    )

    # Apply the wall shade if required AND if I=0
    if wall_shade and not getattr(profile.bfield, "has_i", False):
        _wall_shade(
            canvas=canvas,
            x=theta_lim[0],
            y=getattr(profile, "psi_wall" + suffix).m,
            width=theta_lim[1],
            height=Pthetaspan[1],
        )

    # Energy labels on contour lines (creates gaps to the contour for some
    # reason)
    # plt.clabel(C, inline=1, fontsize=10, zorder=10)

    ax.format_coord = grid_cursor_format_coord(
        x=theta_grid, y=Ptheta_grid, z=W_grid, z_label=E_units
    )

    # Make plot interactive if requested
    if _final:
        plt.ion()
        plt.show(block=True)
        logger.info("\tProfile contour drawn internally.")
    else:
        logger.info("\tProfile contour drawn end shown (end of pipeline) ")


# ---------------------------------------------------------------------------


def _wall_shade(canvas, x, y, width, height):
    """Displays a shade above ψ=ψ_wall."""

    fig, ax = canvas

    rect = Rectangle(
        (x, y),
        width=2 * np.pi,
        height=height,
        alpha=0.2,
        color="k",
    )
    ax.add_patch(rect)
