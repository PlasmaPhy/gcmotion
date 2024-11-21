r"""
Plots the energy contour lines of the :math:`\theta-P_\theta`
plain of the Hamiltonian.

Can also plot the current particle's :math:`\theta-P_\theta` drift.
Should be False when running with multiple initial conditions.

We can optionally remove the :math:`\Phi` term from the Hamiltonian,
to observe the change of the drifts with the added electric field.

The x-axis (angle) limits can be either [-π,π] or [0,2π].

The optional arguements are only used when plotting drifts
from multiple particles in the same canvas.

Example
-------

.. code-block:: python

    gcm.energy_contour(
        cwp, theta_lim = [-np.pi ,np.pi], Ptheta_lim="auto",
        plot_drift=True, contour_Phi=True, units="SI",
        levels=20
    )

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import ticker  # for contour plot locator

from gcmotion.utils.logger_setup import logger
from gcmotion.utils.energy_Ptheta import energy_Ptheta
from gcmotion.utils.plot_utils import yspan

from gcmotion.plotters.drift import drift

from gcmotion.configuration.plot_parameters import energy_contour as config


def energy_contour(
    cwp,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: str | list = "auto",
    plot_drift: bool = True,
    contour_Phi: bool = True,
    units: str = "SI",
    levels: int = None,
    wall_shade: bool = True,
    **params,
):
    r"""
    Plots the energy contour lines of the :math:`\theta-P_\theta`
    plain of the Hamiltonian.

    Can also plot the current particle's :math:`\theta-P_\theta` drift.
    Should be False when running with multiple initial conditions.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    theta_lim : list, optional
        Plot xlim. Must be either [0,2π] or [-π,π]. Defaults to [-π,π].
    psi_lim : list | str, optional
        If a list is passed, it plots between the 2 values of
        :math:`\psi`. If "auto" is passed, it automatically sets
        the optimal :math:`\psi` limits. Defaults to 'auto'.
    plot_drift : bool, optional
        Whether or not to plot :math:`\theta-P_\theta` drift on top.
        Defaults to True.
    contour_Phi : bool, optional
        Whether or not to add the Φ term in the energy contour.
        Defaults to True.
    units : str, optional
        The unit system. Can be either 'SI' or 'NU'. Defauls to "SI".
    levels : int, optional
        The number of contour levels. Defaults to Config setting.
    wall_shade : bool, optional
        Whether to shade the region :math:`\psi/\psi_{wall} > 1`.
        Note that this only makes sense in the case of
        :math:`I=0`, where :math:`P_\theta=\psi`
        Defaults to True.
    params : dict, optional
        Extra plotting parameters:

            * different_colors : bool
                Whether or not not use different colors for every drift.
                Defaults to False.
            * plot_initial: bool
                Whether or not to plot the starting points of each drift.
                Defaults to True.
    """
    if params.get("_internal_call", False):  # dont pop it yet
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    suffix = "NU" if units == "NU" else "" if units == "SI" else ""
    logger.info(f"Plotting contour plot in {
                "NU" if suffix == "NU" else "SI"}...")

    # Get all needed attributes first
    # Use cwp.psiNU to calculate the contour, since psi cannot be solved
    # for Ptheta in the general case if I =/= 0.
    Q = getattr(cwp, "Q")
    psi_wallNU = getattr(cwp, "psi_wallNU").m
    cwppsiNU = getattr(cwp, "psiNU").copy()
    muNU = getattr(cwp, "muNU").magnitude
    PzetaNU = getattr(cwp, "Pzeta0NU").magnitude
    profile = getattr(cwp, "profile")

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!

    # Grab or create canvas
    if canvas is None:
        fig = plt.figure(**config["fig_parameters"])
        ax = fig.add_subplot(111)
        canvas = (fig, ax)
        logger.debug("\tContour: Creating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tContour: Using existing canvas.")

    if theta_lim not in [[0, 2 * np.pi], [-np.pi, np.pi]]:
        print(
            "'theta_lim' must be either [0,2*np.pi] or [-np.pi,np.pi]. Defaulting to [-np.pi,np.pi]")
        logger.warning(
            "\tContour: Invalid 'theta_lim': Defaulting to [-np.pi,np.pi]")
        theta_lim = [-np.pi, np.pi]

    # Plot drift of Ptheta-theta
    if plot_drift:
        drift(
            cwp,
            angle="theta",
            lim=theta_lim,
            units=units,
            _internal_call=True,
            canvas=canvas,
            **params,
        )
        logger.debug("\tContour: Plotting particle's Pθ drift.")

    # Find the psi-grid span of psi upon which to calculate the W and Ptheta
    # grids. Psi must be passed in [NU] despite the "units" arguement.
    if isinstance(psi_lim, str):  # If "auto"
        psi_lim = yspan(cwppsiNU.magnitude, psi_wallNU)
    elif isinstance(psi_lim, (list, np.ndarray)):  # If relative to psi_wall
        psi_lim = np.array(psi_lim) * Q("NUpsi_wall")
        psi_lim = psi_lim.to("NUmagnetic_flux").m

    # Calculate Energy AND Ptheta grid values in [NU] based on the
    # psi limits calculated above.
    grid_density = config["contour_grid_density"]
    theta, psiNU = np.meshgrid(
        np.linspace(theta_lim[0], theta_lim[1], grid_density),
        np.linspace(psi_lim[0], psi_lim[1], grid_density),
    )

    WNU, PthetaNU = energy_Ptheta(
        theta=theta,
        psi=psiNU,
        mu=muNU,
        Pzeta=PzetaNU,
        profile=profile,
        contour_Phi=contour_Phi,
    )

    # Quantify them, convert them to SI if needed and calculate their spans
    WNU = Q(WNU, "NUJoule")
    PthetaNU = Q(PthetaNU, "NUMagnetic_flux")

    if units == "NU":
        Eunits = "NUJoule"
        W_contour = WNU.m
        Ptheta_contour = PthetaNU.m
    elif units == "SI":
        Eunits = "keV"
        W_contour = WNU.to("kev").m
        Ptheta_contour = PthetaNU.to("Magnetic_flux").m

    Wspan = np.array([W_contour.min(), W_contour.max()])
    Pthetaspan = [Ptheta_contour.min(), Ptheta_contour.max()]

    # Configure contour plot
    if levels is None:  # If non is given
        levels = config["contour_levels"]

    if config["locator"] == "log":
        locator = ticker.LogLocator(base=config["log_base"], numticks=levels)
    else:
        locator = ticker.MaxNLocator(nbins=levels)

    contour_kw = {
        "vmin": Wspan[0],
        "vmax": Wspan[1],
        "levels": levels,
        "cmap": config["contour_cmap"],
        "locator": locator,
        "zorder": 1,
    }

    # Contour plot
    # Make sure Ptheta is also in the required units
    C = ax.contourf(theta, Ptheta_contour, W_contour, **contour_kw)
    ax.set_xlabel(r"$\theta$ [rads]")
    ax.set_ylabel(
        rf"$P_\theta$ [{Q(suffix+"Magnetic_flux").units:~P}]", rotation=90
    )
    ticks = ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"]
    plt.xticks(np.linspace(-2 * np.pi, 2 * np.pi, 9), ticks)
    ax.set(xlim=theta_lim, ylim=Pthetaspan)
    ax.set_facecolor("white")

    # Apply the wall shade if required AND if I=0
    if wall_shade and not getattr(cwp.bfield, "has_i", False):
        _wall_shade(
            canvas=canvas,
            x=theta_lim[0],
            y=getattr(cwp, "psi_wall" + suffix).m,
            width=theta_lim[1],
            height=Pthetaspan[1],
        )

    # Energy labels on contour lines (creates gaps to the contour for some
    # reason)
    # plt.clabel(C, inline=1, fontsize=10, zorder=10)

    # Cursor
    Xflat, Yflat, Zflat = (
        theta.flatten(),
        Ptheta_contour.flatten(),
        W_contour.flatten(),
    )

    def fmt(x, y):
        # get closest point with known data
        dist = np.linalg.norm(np.vstack([Xflat - x, Yflat - y]), axis=0)
        idx = np.argmin(dist)
        z = Zflat[idx]
        return f"theta={x:.4g}  Ptheta={y:.5f}  E={z:.4g}[{Eunits}]"

    plt.gca().format_coord = fmt

    # Color bar Plotting
    Esuffix = "NU" if suffix == "NU" else "keV"
    if not _internal_call:  # FIXME
        energies = getattr(cwp, "E" + Esuffix).magnitude
        _ = _cbar(energies=energies, units=units, canvas=canvas, C=C)

    # Make plot interactive if requested
    if not _internal_call:
        fig.set_tight_layout(True)
        plt.ion()
        plt.show(block=True)
        logger.info("--> Energy contour successfully plotted (returned null).")
    else:
        # fmt:skip
        logger.info(
            "--> Energy contour successfully plotted "
            + "(returned contour object)"
        )
        return C


# ---------------------------------------------------------------------------


def _cbar(energies, units, canvas, C):
    """Returns the height of the energy colorbar label.

    Parameters
    ----------
    energies : np.ndarray
        The particle's energies in the correct unit.
    units : str
        The energy units.
    canvas : 2-tuple
        2-tuple containing (fig, ax).
    C : contour object
        The object returned from plt.contour().

    Returns
    -------
        The colorbar object.
    """
    logger.debug("Calling _cbar_energy()")

    fig, ax = canvas
    energies = list([energies])

    cbar = fig.colorbar(
        C,
        ax=ax,
        fraction=0.03,
        pad=0.1,
        label=f"E[{units}]",
        format="{x:.4g}",
        extendrect="auto",
    )
    cbar_kw = {
        "linestyle": "-",
        "zorder": 3,
        "color": config["cbar_color"],
    }
    for E in energies:
        cbar.ax.plot([0, 1], [E, E], **cbar_kw)

    return cbar


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
