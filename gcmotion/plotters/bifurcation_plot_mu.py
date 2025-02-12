r"""
Calculates and plots the bifurcation diagram of the fixed points for multiple profiles
with different :math:`\mu`'s.
"""

import matplotlib.pyplot as plt
import numpy as np
from time import time
from gcmotion.utils.logger_setup import logger

from gcmotion.scripts.bifurcation import bifurcation
from collections import deque
from gcmotion.configuration.fixed_points_bifurcation_parameters import BifurcationPlotConfig


def bifurcation_plot_mu(profiles: list | deque, **kwargs):
    r"""Draws the bifurcation diagrams for the :math:`theta`'s  fixed,
    the :math:`P_{theta}`'s fixed and the number of fixed points found for
    each :math:`P_{\zeta}`.

    :meta public:

        Parameters
        ----------
        profiles : list, deque
            List of profile objects that contain Tokamak and Particle information.
        Other Parameters
        ----------
        thetalim : list, optional
            Limits of the of the :math:`\theta`, :math:`\psi` search area with respect
            to the :math:`\theta` variable. Defaults to [-:math:`\pi`, :math:`\pi`].
        psilim : list, optional
            Limits of the of the :math:`\theta`, :math:`\psi` search area with respect
            to the :math:`\psi` variable. Defaults to [0.01 , 1.8]. CUTION: The limits are given
            normalized to :math:`\psi_{wall}`.
        method : str, optional
            String that indicates which method will be used to find the systems fixed
            points in :py:func:`single_fixed_point`. Can either be "fsolve" (deterministic)
            or "differential evolution" (stochastic). Defaults to "fsolve".
        dist_tol : float, optional
            Tolerance below which two fixed points are not considered distinct. The differences between
            both :math:`\theta` and :math:`\psi` of the fixed points must be below this tolerance for
            the fixed points to be considered the same. Defaults to 1e-3.
        fp_ic_scan_tol : float, optional
            Tolerance below which the sum of the squares of the time derivatives of the
            :math:`\theta` and :math:`\psi` variavles is considered zero. It is passed into
            :py:func:`fp_ic_scan`. Defaults to 5 * 1e-8.
        ic_theta_grid_density : int, optional
            Density of the :math:`\theta`, :math:`\psi` 2D grid to be scanned for initial conditiond
            (fixed points candidates) with respect to the :math:`\theta` variable. It is passed into
            :py:func:`fp_ic_scan` Defaults to 400.
        ic_psi_grid_density : int, optional
            Density of the :math:`\theta`, :math:`\psi` 2D grid to be scanned for initial conditiond
            (fixed points candidates) with respect to the :math:`\psi` variable. It is passed into
            :py:func:`fp_ic_scan` Defaults to 400.
        random_fp_init_cond : bool, optional
            Boolean determining weather random initial conditions are to be used instead of those
            provided by :py:func:`fp_ic_scan`. Defaults to ``False``.
        fp_info : bool, optional
            Boolean determining weather fixed points' information is to be printed. Defaults to ``False``.
        bif_info: bool, optional
            Boolean that determines weather information regarding the bifurcation process is to
            be printed. Defaults to ``False``.
        fp_ic_info : bool, optional
            Boolean determing weather information on the initial condition is to be printed.
            Defaults to ``False``.
        plot_energy_bif : bool, optional
            Boolean determining weather the energy of each fixed point of each profile (each :math:`\P_{\zeta}`)
            is to be plotted. Defaults to ``False``.
        energy_units : str, optional
            String specifying the unit of the calculated fixed points' energies. Defaults to ``"NUJoule"``.
        energies_info : bool, optional
            Boolean determining weather information on the fixed points' energies is to be printed.
            Defaults to ``True``.
        fp_LAR_thetas : bool, optional
            Boolean determining weather the theta values for which fixed points occur are to be
            considered known (LAR thetas are 0 and :math:`\pi`). Defaults to ``False``.
        fp_only_confined : bool, optional
            Boolean determining if the search for :math:`\psi_{fixed}` will be conducted only for
            :math:`\psi` < :math:`\psi_{wall}` (confined particles). Defaults to ``False``.
    """

    # Unpack parameters
    config = BifurcationPlotConfig()
    for key, value in kwargs.items():
        setattr(config, key, value)

    start = time()
    # CAUTION: The bifurcation function takes in psis_fixed but returns P_thetas_fixed
    X_thetas, X_P_thetas, O_thetas, O_P_thetas, num_of_XP, num_of_OP, X_energies, O_energies = (
        bifurcation(
            profiles=profiles,
            calc_energies=config.plot_energy_bif,
            which_COM="mu",
            **kwargs,
        )
    )

    print(f"BIFURCATION RUN IN {(time() - start)/60:.1f} mins")

    profile1 = profiles[0]
    profileN = profiles[-1]
    logger.info(
        f"Ran bifurcation script for bifurcation plot with N={len(profiles)} mus={profile1.muNU}...{profileN.muNU}"
    )
    psi_wallNU = profile1.psi_wall.to("NUMagnetic_flux")
    P_theta_wallNU = profile1.findPtheta(psi=psi_wallNU, units="NUCanonical_momentum").m

    fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    plt.xlabel(r"${\mu}$" + f"[{profiles[1].muNU.units}]")
    fig.suptitle("Fixed Points Bifurcation Diagram")

    ax_theta = ax[0]
    ax_P_theta = ax[1]
    ax_ndfp = ax[2]

    mu_plot = []

    X_theta_plot = []
    X_P_theta_plot = []

    O_theta_plot = []
    O_P_theta_plot = []

    # Theta Fixed Bifurcation
    for i, profile in enumerate(profiles):
        mu = profile.muNU.m
        y_list = X_thetas[i]
        mu_plot.extend([mu] * len(list(y_list)))
        X_theta_plot.extend(y_list)

    ax_theta.set_ylabel(r"$\theta_s$ Fixed")
    ax_theta.set_yticks([-np.pi, 0, np.pi])
    ax_theta.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax_theta.set_ylim([-np.pi - 0.5, np.pi + 0.5])
    ax_theta.scatter(mu_plot, X_theta_plot, s=2, color="#E65100")

    logger.info(f"Made Xthetas fixed bifurcation plot")

    mu_plot = []

    for i, profile in enumerate(profiles):
        mu = profile.muNU.m
        y_list = O_thetas[i]
        mu_plot.extend([mu] * len(list(y_list)))
        O_theta_plot.extend(y_list)

    ax_theta.scatter(mu_plot, O_theta_plot, s=2)

    logger.info(f"Made Othetas fixed bifurcation plot")

    mu_plot1 = []

    # P_theta Fixed Bifurcation
    for i, profile in enumerate(profiles):
        mu = profile.muNU.m
        y_list = X_P_thetas[i]
        mu_plot1.extend([mu] * len(list(y_list)))
        X_P_theta_plot.extend(y_list)

    mu_plot2 = []

    for i, profile in enumerate(profiles):
        mu = profile.muNU.m
        y_list = O_P_thetas[i]
        mu_plot2.extend([mu] * len(list(y_list)))
        O_P_theta_plot.extend(y_list)

    # Set the upper limit of the y axis properly
    ul = 1.05 * P_theta_wallNU

    # Combine the two lists for comparison
    combined_P_theta_plot = X_P_theta_plot + O_P_theta_plot

    if max(combined_P_theta_plot) > P_theta_wallNU:
        ul = max(combined_P_theta_plot) * 1.05

    ax_P_theta.set_ylabel(r"$P_{\theta_s}$ Fixed")
    ax_P_theta.set_ylim(0, ul)
    ax_P_theta.scatter(mu_plot1, X_P_theta_plot, s=2, color="#E65100", label="X points")
    ax_P_theta.scatter(mu_plot2, O_P_theta_plot, s=2, label="O points")
    ax_P_theta.axhline(y=P_theta_wallNU, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax_P_theta.legend(loc="lower left")

    logger.info(f"Made P_thetas fixed bifurcation plot")

    # Number of distinct fixed points Diagram
    mus = [profile.muNU.m for profile in profiles]
    ax_ndfp.set_ylabel("Number of Fixed Points")
    ax_ndfp.scatter(mus, num_of_XP, s=2, color="#E65100", label="X points")
    ax_ndfp.scatter(mus, num_of_OP, s=2, label="O points")
    ax_ndfp.legend()

    logger.info(f"Made number of fixed points bifurcation plot")

    if config.plot_energy_bif:
        fig, ax = plt.subplots(1, 1, figsize=(9, 7), sharex=True)
        plt.xlabel(r"$E-{\mu}B_0$" + f"[{config.energy_units}]")
        ax.set_ylabel(r"$\mu$" + f"[{profile1.muNU.units}]")

        X_energies_plot = []
        O_energies_plot = []

        mu_plot1 = []

        # XPoints Energies - muB0 Bifurcation
        for i, profile in enumerate(profiles):
            mu = profile.muNU
            muB0 = (mu * profile.Q(1, "NUTesla")).to(config.energy_units)
            y_list = X_energies[i].m - muB0.m

            # Ensure y_list is iterable
            if np.isscalar(y_list):
                y_list = [y_list]

            mu_plot1.extend([mu.m] * len(list(y_list)))
            X_energies_plot.extend(y_list)

        mu_plot2 = []

        # OPoints Energies - muB0 Bifurcation
        for i, profile in enumerate(profiles):
            mu = profile.muNU
            muB0 = (mu * profile.Q(1, "NUTesla")).to(config.energy_units)
            y_list = O_energies[i].m - muB0.m

            # Ensure y_list is iterable
            if np.isscalar(y_list):
                y_list = [y_list]

            mu_plot2.extend([mu.m] * len(list(y_list)))
            O_energies_plot.extend(y_list)

        ax.scatter(X_energies_plot, mu_plot1, s=2, color="#E65100", label="X points")
        ax.scatter(O_energies_plot, mu_plot2, s=2, label="O points")

        logger.info(f"Made fixed points' energies bifurcation plot")

        ax.legend()

    plt.ion()
    plt.show(block=True)
