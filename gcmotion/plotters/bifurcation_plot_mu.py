r"""
Calculates and plots the bifurcation diagram of the fixed points for multiple profiles
with different :math:`\mu`'s.
"""

import matplotlib.pyplot as plt
from time import time
from gcmotion.utils.logger_setup import logger

from gcmotion.plotters._bif_base._thetas_bif_plot import _thetas_bif_plot
from gcmotion.plotters._bif_base._P_thetas_bif_plot import _P_thetas_bif_plot
from gcmotion.plotters._bif_base._ndfp_bif_plot import _ndfp_bif_plot
from gcmotion.plotters._bif_base._tpb_plot import _plot_trapped_passing_boundary

from gcmotion.scripts.bifurcation import bifurcation
from collections import deque
from gcmotion.configuration.plot_parameters import BifurcationPlotConfig


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

    which_COM_loc = "mu"

    start = time()
    # CAUTION: The bifurcation function takes in psis_fixed but returns P_thetas_fixed
    X_thetas, X_P_thetas, O_thetas, O_P_thetas, num_of_XP, num_of_OP, X_energies, O_energies = (
        bifurcation(
            profiles=profiles,
            calc_energies=config.plot_energy_bif,
            which_COM=which_COM_loc,
            **kwargs,
        )
    )

    print(f"BIFURCATION RUN IN {(time() - start)/60:.1f} mins")

    profile1 = profiles[0]
    profileN = profiles[-1]
    logger.info(
        f"Ran bifurcation script for bifurcation plot with N={len(profiles)}, Pz={profile1.PzetaNU} mus={profile1.muNU}...{profileN.muNU}"
    )

    # Create figure
    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
        "facecolor": config.facecolor,
        "sharex": config.sharex,
    }

    fig, ax = plt.subplots(3, 1, **fig_kw)
    plt.xlabel(r"${\mu}$" + f"[{profile1.muNU.units}]")
    fig.suptitle("Fixed Points Bifurcation Diagram")

    ax_theta = ax[0]
    ax_P_theta = ax[1]
    ax_ndfp = ax[2]

    ax_theta.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax_P_theta.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax_ndfp.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    # Fixed thetas bifurcation diagram
    _thetas_bif_plot(
        profiles=profiles,
        X_thetas=X_thetas,
        O_thetas=O_thetas,
        which_COM=which_COM_loc,
        ax=ax_theta,
    )

    logger.info(f"Made Xthetas, Othetas fixed bifurcation plot")

    # P_theta Fixed Bifurcation
    _P_thetas_bif_plot(
        profiles=profiles,
        X_P_thetas=X_P_thetas,
        O_P_thetas=O_P_thetas,
        which_COM=which_COM_loc,
        ax=ax_P_theta,
    )

    logger.info(f"Made P_thetas fixed bifurcation plot")

    # Number of distinct fixed points Diagram
    _ndfp_bif_plot(
        profiles=profiles,
        num_of_XP=num_of_XP,
        num_of_OP=num_of_OP,
        which_COM=which_COM_loc,
        ax=ax_ndfp,
    )

    logger.info(f"Made number of fixed points bifurcation plot")

    if config.plot_energy_bif:
        _plot_trapped_passing_boundary(
            profiles=profiles,
            X_energies=X_energies,
            O_energies=O_energies,
            which_COM=which_COM_loc,
            config=config,
            ax=ax,
        )

        logger.info(f"Made fixed points' energies bifurcation plot")

    plt.ion()
    plt.show(block=True)
