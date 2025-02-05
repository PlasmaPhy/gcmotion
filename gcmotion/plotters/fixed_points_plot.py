r"""
Draws fixed points plot. This method is called internally by ``profile_contour()``
as well.
"""

from gcmotion.scripts.fixed_points import fixed_points as fp
from gcmotion.entities.profile import Profile

from gcmotion.utils.XO_points_classification import XO_points_classification as xoc
from gcmotion.utils.points_psi_to_P_theta import points_psi_to_P_theta

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np

from time import time

from gcmotion.utils.logger_setup import logger


def fixed_points_plot(
    profile: Profile,
    theta_lim: list,
    psi_lim: list,
    method: str = "fsolve",
    dist_tol: float = 1e-3,
    fp_ic_scan_tol: float = 5 * 1e-8,
    ic_theta_grid_density: int = 400,
    ic_psi_grid_density: int = 400,
    random_init_cond: bool = False,
    info: bool = False,
    ic_info: bool = False,
    plot_init_cond: bool = False,
    LAR_thetas: bool = False,
    only_confined: bool = False,
    ax: Axes = None,
    **args,
):
    r"""

        :meta public:

    Function that creates a plot of the GC Hamiltonian's fixed points. Most of its arguments
    will be passed into :py:func:`fixed_points`. This function is almost always used in
    :py:func:`profile_contour`, in order to visualize the results.

    Parameters
            ----------
            profile : Profile
                Profile object that contains Tokamak and Particle information.
            theta_lim : list, optional
                Limits of the of the :math:`\theta`, :math:`\psi` search area with respect
                to the :math:`\theta` variable. Defaults to [-:math:`\pi`, :math:`\pi`].
            method : str, optional
                String that indicates which method will be used to find the systems fixed
                points in :py:func:`single_fixed_point`. Can either be "fsolve" (deterministic)
                or "differential evolution" (stochastic). Defaults to "fsolve".
            psi_lim : list, optional
                Limits of the of the :math:`\theta`, :math:`\psi` search area with respect
                to the :math:`\psi` variable. Defaults to [0.01 , 1.8]. CUTION: The limits are given
                normalized to :math:`\psi_{wall}`.
            dist_tol : float
                Tolerance below which two fixed points are not considered distinct. The differences between
                both :math:`\theta` and :math:`\psi` of the fixed points must be below this tolerance for
                the fixed points to be considered the same. Defaults to 1e-3.
            fp_ic_scan_tol : float
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
            random_init_cond : bool, optional
                Boolean determining weather random initial conditions are to be used instead of those
                provided by :py:func:`fp_ic_scan`. Defaults to ``False``.
            info : bool, optional
                Boolean determining weather fixed points' information is to be printed. Defaults to ``False``.
            ic_info : bool, optional
                Boolean determing weather information on the initial condition is to be printed.
                Defaults to ``False``.
            plot_init_cond : bool, optional
                Boolean that determines weather the initial conditions passed into :py:func:`fixed_points`
                will be plotted. Defaults to ``False``.
            LAR_thetas : bool, optional
                Boolean determining weather the theta values for which fixed points occur are to be
                considered known (LAR thetas are 0 and :math:`\pi`). Defaults to ``False``.
            only_confined : bool, optional
                Boolean determining if the search for :math:`\psi_{fixed}` will be conducted only for
                :math:`\psi` < :math:`\psi_{wall}` (confined particles). Defaults to ``False``.
            ax : Axes, optional
                The axes upon which the plot is to be realized.
            args : dict, optional
                Extra arguements.
    """

    # Unpack params
    _internal_call = args.pop("_internal_call", False)  # POP!

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    logger.info(f"Plotting fixed points")

    start = time()
    # Calculate fixed points
    _, fixed_points, initial_conditions = fp(
        method=method,
        profile=profile,
        Q=profile.Q,
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        dist_tol=dist_tol,
        fp_ic_scan_tol=fp_ic_scan_tol,
        ic_theta_grid_density=ic_theta_grid_density,
        ic_psi_grid_density=ic_psi_grid_density,
        random_init_cond=random_init_cond,
        info=info,
        ic_info=ic_info,
        LAR_thetas=LAR_thetas,
        only_confined=only_confined,
    )
    print(f"\n FIXED POINTS RUN IN {(time() - start):.1f}s\n")

    # CAUTION: The xoc function takes in psis_fixed but returns P_thetas_fixed
    X_points, O_points = xoc(unclassified_fixed_points=fixed_points, profile=profile)

    # Check axes
    if ax is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        logger.debug("\tCreating a new canvas.")

    # Convert deque to numpy arrays for easy manipulation
    X_thetas, X_P_thetasNU = zip(*X_points) if X_points else ([], [])
    O_thetas, O_P_thetasNU = zip(*O_points) if O_points else ([], [])

    X_P_thetas = profile.Q(X_P_thetasNU, "NUmagnetic_flux").to("Magnetic_flux")
    O_P_thetas = profile.Q(O_P_thetasNU, "NUmagnetic_flux").to("Magnetic_flux")

    ax.set_xticks([-np.pi, 0, np.pi])
    ax.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax.set_xlim([-np.pi, np.pi])

    if len(X_thetas) > 0 and len(X_P_thetas) > 0:
        ax.scatter(X_thetas, X_P_thetas, marker="x", color="#80FF80", s=100)

    if len(O_thetas) > 0 and len(O_P_thetas) > 0:
        ax.scatter(O_thetas, O_P_thetas, marker="o", edgecolor="yellow", facecolors="none", s=100)

    if plot_init_cond:

        # Turn points (theta, psi) --> (theta, P_theta)
        initial_conditions = points_psi_to_P_theta(initial_conditions, profile=profile)

        thetas_init, P_thetas_initNU = zip(*initial_conditions) if initial_conditions else ([], [])

        P_thetas_init = profile.Q(P_thetas_initNU, "NUmagnetic_flux").to("Magnetic_flux")
        ax.scatter(thetas_init, P_thetas_init.m, marker=">", color="red", s=100)
