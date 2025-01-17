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
    method: str = "differential evolution",
    dist_tol: float = 1e-3,
    fp_ic_scan_tol: float = 5 * 1e-8,
    ic_theta_grid_density: int = 1000,
    ic_psi_grid_density: int = 1000,
    random_init_cond: bool = False,
    info: bool = False,
    ic_info: bool = False,
    plot_init_cond: bool = False,
    LAR_thetas: bool = False,
    ax: Axes = None,
    **args,
):
    r"""Draws fixed points plot.

    This method is called internally by ``energy_contour()``
    as well.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    theta_lim : list
        List containing the limits for :math:`\theta` (x- axis limits).
    psi_lim : list
        List containing the limits for :math:`\psi` and
        consequently :math:`\psi` (y- axis limits).
    dist_tol : float, optional
        Tolerance that determines distinct fixed points. If both :math:`\theta` and
        :math:`\psi` elements of a fixed point are less than :py:data:`dist_tol` apart
        the two fixed points are not considered distinct.
    ic_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\theta` variable
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution`
        will be conducted. Will be passed to :py:func:`fixed_points`.
    ic_psi_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\psi` variable
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution`
        will be conducted.  Will be passed to :py:func:`fixed_points`.
    info : bool, optional
        Passed into ``fixed_poits()``. Determines weather to print information
        about the fixed points (number, values). Defaults to ``False``.
    params : dict, optional
        Extra arguements if called for many particles:
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
    )
    print(f"\n FIXED POINTS RUN IN {(time() - start):.1f}s\n")

    # CAUTION: The xoc function takes in psis_fixed but returns P_thetas_fixed
    X_points, O_points = xoc(
        unclassified_fixed_points=fixed_points,
        profile=profile,
    )

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
