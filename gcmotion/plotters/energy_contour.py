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
        cwp, lim = [-np.pi ,np.pi], psi_lim="auto", 
        plot_drift=True, contour_Phi=True, units="keV", 
        levels=20
    )

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import ticker  # for contour plot locator

from gcmotion.utils._logger_setup import logger

from gcmotion.plotters.drift import drift

from gcmotion.configuration.plot_parameters import energy_contour as config


def energy_contour(
    cwp,
    lim: list = [-np.pi, np.pi],
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
    lim : list, optional
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
        The unit system. Can be either 'SI' or 'NU'. Defauls to "SI".
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
    """
    suffix = "NU" if units == "NU" else "" if units == "SI" else ""
    logger.info(f"Plotting contour plot in {"NU" if suffix=="NU" else "SI"}...")  # fmt: skip

    # Get all needed attributes first
    psi_wall = getattr(cwp, "psi_wallNU").copy()
    cwp_psi = getattr(cwp, "psiNU").copy()
    Q = getattr(cwp, "Q")

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    # Grab or create canvas
    if canvas is None:
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111)
        canvas = (fig, ax)
        logger.debug("\tContour: Creating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tContour: Using existing canvas.")

    # Plot drift if requested
    if plot_drift:
        drift(
            cwp,
            angle="theta",
            lim=lim,
            _internal_call=True,
            canvas=canvas,
            **params,
        )
        logger.debug("\tContour: Plotting particle's Pθ drift.")

    # Set psi limits (Normalised to psi_wall)
    if type(psi_lim) is str:
        if psi_lim == "auto":
            psi_diff = cwp_psi.max() - cwp_psi.min()
            psi_mid = (cwp_psi.max() + cwp_psi.min()) / 2
            psi_diff, psi_mid = psi_diff.magnitude, psi_mid.magnitude
            psi_lower = max(0, psi_mid - 0.6 * psi_diff)
            psi_higher = psi_mid + 0.6 * psi_diff
            psi_lim = np.array([psi_lower, psi_higher])
    else:
        psi_lim = np.array(psi_lim) * psi_wall

    psi_lim = Q(psi_lim, "NUMagnetic_flux")
    psi_min = psi_lim[0].magnitude
    psi_max = psi_lim[1].magnitude

    # Set theta lim.
    theta_min, theta_max = lim

    # Calculate Energy values
    grid_density = config["contour_grid_density"]
    theta, psi = np.meshgrid(
        np.linspace(theta_min, theta_max, grid_density),
        np.linspace(psi_min, psi_max, grid_density),
    )
    WNU = Q(_calcWNU_grid(cwp, theta, psi, contour_Phi, units),"NUJoule")  # fmt: skip
    W = (
        WNU.to("keV")
        if units == "SI"
        else WNU if units == "NU" else WNU.to("keV")
    )
    span = np.array([W.magnitude.min(), W.magnitude.max()])

    # Configure contour plot
    if levels is None:  # If non is given
        levels = config["contour_levels"]

    if config["locator"] == "log":
        locator = ticker.LogLocator(base=1.05, numticks=levels)
    else:
        locator = ticker.MaxNLocator(nbins=levels)

    contour_kw = {
        "vmin": span[0],
        "vmax": span[1],
        "levels": levels,
        "cmap": config["contour_cmap"],
        "locator": locator,
        "zorder": 1,
    }

    # Contour plot
    C = ax.contourf(theta, psi / psi_wall.magnitude, W.magnitude, **contour_kw)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\psi/\psi_{wall}$", rotation=90)
    ticks = ["-2π", "-3π/2", "-π", "-π/2", "0", "π/2", "π", "3π/2", "2π"]
    plt.xticks(np.linspace(-2 * np.pi, 2 * np.pi, 9), ticks)
    ax.set(
        xlim=[theta_min, theta_max],
        ylim=np.array(psi_lim.magnitude) / psi_wall.magnitude,
    )
    ax.set_facecolor("white")

    # Wall shade
    if wall_shade:  # ψ_wall boundary rectangle
        rect = Rectangle(
            (lim[0], 1),
            2 * np.pi,
            psi_max / psi_wall.magnitude,
            alpha=0.2,
            color="k",
        )
        ax.add_patch(rect)

    # Energy labels on contour lines (creates gaps to the contour for some reason)
    # plt.clabel(C, inline=1, fontsize=10, zorder=10)

    # Cursor
    def fmt(theta, psi):
        E = _calcWNU_grid(
            cwp,
            theta,
            psi * cwp.psi_wall,
            cwp.Pzeta0,
            contour_Phi=contour_Phi,
            units=units,
        )
        return (
            "theta={theta:.4f},\tpsi={psi:.4f},\tE={E:.4f} ".format(
                theta=theta, psi=psi, E=E
            )
            + units
        )

    plt.gca().format_coord = fmt

    # Color bar Plotting
    if not _internal_call:
        energies = cwp.__getattribute__("E" + suffix).magnitude
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


def _calcWNU_grid(
    cwp,
    theta: np.array,
    psi: np.array,
    contour_Phi: bool,
    units: str,
):
    r"""Returns a single value or a grid of the calculated Hamiltonian.

    Only to be called internally, by ``energy_contour()``..

    Args:
        theta (np.array): The :math:`\theta` values.
        psi (np.array): The :math:`\psi` values.
        contour_Phi (bool): Whether or not to add the electric potential term
            :math:`e\Phi`.
        units (str): The energy units.
    """
    # Get all needed attributes first
    mu = cwp.muNU.magnitude
    Pzeta = cwp.Pzeta0NU.magnitude
    qfactor = cwp.qfactor
    bfield = cwp.bfield
    efield = cwp.efield
    g = bfield.gNU.magnitude

    r = np.sqrt(2 * psi)
    B = bfield.b(r, theta)
    psip = qfactor.psip_of_psi(psi)

    rho = (Pzeta + psip) / g**2
    W = (1 / 2) * rho**2 * B**2 + mu * B  # Without Φ

    # Add Φ if asked
    if contour_Phi:
        Phi = efield.Phi_of_psi(psi)
        W += Phi  # all normalized
        logger.debug("\tContour: Adding Φ term to the contour.")

    return W


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
