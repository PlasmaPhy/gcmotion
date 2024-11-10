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
    Ptheta_lim: str | list = "auto",
    plot_drift: bool = True,
    contour_Phi: bool = True,
    units: str = "SI",
    levels: int = None,
    wall_shade: bool = False,
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
    Ptheta_lim : list | str, optional
        If a list is passed, it plots between the 2 values of
        :math:`P_\theta`. If "auto" is passed, it automatically sets
        the optimal :math:`P_\theta` limits. Defaults to 'auto'.
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
        Defaults to False.
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
    logger.info(f"Plotting contour plot in {"NU" if suffix=="NU" else "SI"}...")  # fmt: skip

    # Get all needed attributes first
    # Use cwp.psiNU to calculate the contour, since psi cannot be solved
    # for Ptheta in the general case if I =/= 0.
    Q = getattr(cwp, "Q")
    cwpPtheta = getattr(cwp, "Ptheta" + suffix).copy()
    cwppsiNU = getattr(cwp, "psiNU").copy()  # NU always
    muNU = getattr(cwp, "muNU").magnitude
    PzetaNU = getattr(cwp, "Pzeta0NU").magnitude
    profile = getattr(cwp, "profile")

    # Hard y lim, in the case that the particle escapes
    bound = config["hardylim"] * getattr(cwp, "psi_wall" + suffix)

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

    # Find the psi-grid span of psi upon which to calculate the energy.
    # Psi is in [NU] despite the "units" arguement.
    if isinstance(Ptheta_lim, list):
        pass  # FIXME: find a way to set Ptheta limits for the contour plot
    elif isinstance(Ptheta_lim, str):
        psiNUspan, _ = yspan(cwppsiNU.magnitude)

    # Calculate Energy AND Ptheta grid values in [NU] based on the
    # psi limits calculated above.
    psiNU_min, psiNU_max = (psiNUspan[0], psiNUspan[1])

    grid_density = config["contour_grid_density"]
    theta, psiNU = np.meshgrid(
        np.linspace(theta_lim[0], theta_lim[1], grid_density),
        np.linspace(psiNU_min, psiNU_max, grid_density),
    )

    WNU, PthetaNU = energy_Ptheta(
        theta=theta,
        psi=psiNU,
        mu=muNU,
        Pzeta=PzetaNU,
        profile=profile,
        contour_Phi=contour_Phi,
    )
    # Convert them to SI if needed and calculate their spans
    WNU = Q(WNU, "NUJoule")
    PthetaNU = Q(PthetaNU, "NUMagnetic_flux")
    if units == "NU":
        W_contour = (
            WNU.m if units == "NU" else WNU.to("keV").m
        )  # Purely numeric
        Ptheta_contour = PthetaNU.m
    elif units == "SI":
        W_contour = WNU.to("kev").m
        Ptheta_contour = PthetaNU.to("Magnetic_flux").magnitude

    Wspan = np.array([W_contour.min(), W_contour.max()])
    Pthetaspan = (
        [Ptheta_contour.min(), Ptheta_contour.max()]
        if Ptheta_lim == "auto"
        else Ptheta_lim
    )

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
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$P_\theta$", rotation=90)
    ticks = ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"]
    plt.xticks(np.linspace(-2 * np.pi, 2 * np.pi, 9), ticks)
    ax.set(xlim=theta_lim, ylim=Pthetaspan)
    ax.set_facecolor("white")

    # Wall shade
    if wall_shade:  # ψ_wall boundary rectangle
        _wall_shade(
            canvas=canvas,
            x=theta_lim[0],
            y=getattr(cwp, "psi_wall" + suffix).m,
            height=Pthetaspan[1],
        )

    # Energy labels on contour lines (creates gaps to the contour for some reason)
    # plt.clabel(C, inline=1, fontsize=10, zorder=10)

    # # Cursor
    # if units == "NU":
    #     cursor_units = "NUJoule"
    #     cursor_Ptheta_units = "NUMagnetic_flux"
    # elif units == "SI":
    #     cursor_units = "keV"
    #     cursor_Ptheta_units = "Magnetic_flux"

    # def fmt(theta, Ptheta):  # FIXME
    #     Ptheta = Q(Ptheta, cursor_Ptheta_units)
    #     E = Q(
    #         energy_Ptheta(
    #             theta=theta,
    #             psi=psiNU,
    #             mu=muNU,
    #             Pzeta=PzetaNU,
    #             profile=profile,
    #             contour_Phi=contour_Phi,
    #         ),
    #         "NUJoule",
    #     ).to(cursor_units)
    #     return "theta={theta:.4g},\tPtheta={Ptheta:.4g~P},\tE={E:.4g~P} \n".format(
    #         theta=theta, Ptheta=Ptheta, E=E
    #     )

    # plt.gca().format_coord = fmt

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
        logger.info("--> Energy contour successfully plotted (returned contour object)\n")  # fmt:skip
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


def _wall_shade(canvas, x, y, height):
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
