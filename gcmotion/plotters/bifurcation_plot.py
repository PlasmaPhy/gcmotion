r"""
Calculates and plots the bifurcation diagram of the fixed points for a Collection 
of particles.

Example
-------

This is how :py:func:`bifurcation_plot` can be called inside a .py script or a .ipynb notebook:

.. code-block:: python

    bifurcation_plot(
        collection=collection,
        theta_lim=[-1.01 * np.pi, 1.01 * np.pi],
        psi_lim=[0.01, 1.3],
        iterations=30,
        info=True,
    )


"""

import matplotlib.pyplot as plt
import numpy as np
from time import time

from gcmotion.scripts.bifurcation import bifurcation
from gcmotion.classes.collection import Collection


def bifurcation_plot(
    collection: Collection,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: list = [0.01, 1.3],
    fp_method: str = "differential evolution",
    dist_tol: float = 1e-3,
    fp_ic_scan_tol: float = 5 * 1e-8,
    ic_theta_grid_density: int = 1000,
    ic_psi_grid_density: int = 1000,
    random_fp_init_cond: bool = False,
    fp_info: bool = False,
    bif_info: bool = False,
    fp_ic_info: bool = False,
):
    r"""Draws the bifurcation diagrams for the :math:`theta`'s  fixed,
    the :math:`P_{theta}`'s fixed and the number of fixed points found for
    each :math:`P_{\zeta0}`.

    :meta public:

    Parameters
    ----------
    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles

    theta_lim : list, optional
        Provides the limits for the solution search area for fixed points
        with regards to the :math:`\theta` variable. It will be passed into
        :py:func:`bifurcation`
    psi_lim : list, optional
        Provides the limits (divided by psi_wall) for the solution search area with regards
        to the :math:`\psi` variable. It will be passed into the "bounds" argument of
        :py:func:`bifurcation`.
    dist_tol : float, optional
        Tolerance that determines distinct fixed points. If both :math:`\theta` and
        :math:`\psi` elements of a fixed point are less than :py:data:`dist_tol` apart
        the two fixed points are not considered distinct.
    ic_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\theta` variable
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution`
        will be conducted. Will be passed to :py:func:`bifurcation`.
    ic_psi_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\psi` variable
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution`
        will be conducted.  Will be passed to :py:func:`bifurcation`.
    info : bool, optional
        Boolean that dictates weather the :math:`P_{\zeta0}` of the particle whose
        fixed points have just been calculated, will be printed alongside the fixed points
        found.
    """

    start = time()
    # CAUTION: The bifurcation function takes in psis_fixed but returns P_thetas_fixed
    X_thetas, X_P_thetas, O_thetas, O_P_thetas, num_of_XP, num_of_OP = bifurcation(
        collection=collection,
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        method=fp_method,
        dist_tol=dist_tol,
        fp_ic_scan_tol=fp_ic_scan_tol,
        ic_theta_grid_density=ic_theta_grid_density,
        ic_psi_grid_density=ic_psi_grid_density,
        random_fp_init_cond=random_fp_init_cond,
        fp_info=fp_info,
        bif_info=bif_info,
        fp_ic_info=fp_ic_info,
    )

    print(f"BIFURCATION RUN IN {(time() - start)/60:.1f} mins")

    particles = collection.particles
    p1 = particles[0]
    psi_wallNU = p1.psi_wallNU.magnitude

    fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    plt.xlabel(r"$P_{\zeta}$ [NUmf]")
    fig.suptitle("Fixed Points Bifurcation Diagram")

    ax_theta = ax[0]
    ax_P_theta = ax[1]
    ax_ndfp = ax[2]

    P_zeta_plot = []

    X_theta_plot = []
    X_P_theta_plot = []

    O_theta_plot = []
    O_P_theta_plot = []

    particles = collection.particles

    # Theta Fixed Bifurcation
    for i, p in enumerate(particles):
        P_zeta = p.Pzeta0NU
        y_list = X_thetas[i]
        P_zeta_plot.extend([P_zeta] * len(list(y_list)))
        X_theta_plot.extend(y_list)

    ax_theta.set_ylabel(r"$\theta_s$ Fixed")
    ax_theta.set_yticks([-np.pi, 0, np.pi])
    ax_theta.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax_theta.set_ylim([-np.pi - 0.5, np.pi + 0.5])
    ax_theta.scatter(P_zeta_plot, X_theta_plot, s=2, color="#E65100")

    P_zeta_plot = []

    for i, p in enumerate(particles):
        P_zeta = p.Pzeta0NU
        y_list = O_thetas[i]
        P_zeta_plot.extend([P_zeta] * len(list(y_list)))
        O_theta_plot.extend(y_list)

    ax_theta.scatter(P_zeta_plot, O_theta_plot, s=2)

    P_zeta_plot1 = []

    # P_theta Fixed Bifurcation
    for i, p in enumerate(particles):
        P_zeta = p.Pzeta0NU
        y_list = X_P_thetas[i]
        P_zeta_plot1.extend([P_zeta] * len(list(y_list)))
        X_P_theta_plot.extend(y_list)

    P_zeta_plot2 = []

    # P_theta Fixed Bifurcation
    for i, p in enumerate(particles):
        P_zeta = p.Pzeta0NU
        y_list = O_P_thetas[i]
        P_zeta_plot2.extend([P_zeta] * len(list(y_list)))
        O_P_theta_plot.extend(y_list)

    # Set the upper limit of the y axis properly
    ul = 1.05 * psi_wallNU

    # Combine the two lists for comparison
    combined_P_theta_plot = X_P_theta_plot + O_P_theta_plot

    if max(combined_P_theta_plot) > psi_wallNU:
        ul = max(combined_P_theta_plot) * 1.05

    ax_P_theta.set_ylabel(r"$P_{\theta_s}$ Fixed")
    ax_P_theta.set_ylim(0, ul)
    ax_P_theta.scatter(P_zeta_plot1, X_P_theta_plot, s=2, color="#E65100", label="X points")
    ax_P_theta.scatter(P_zeta_plot2, O_P_theta_plot, s=2, label="O points")
    ax_P_theta.axhline(y=psi_wallNU, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax_P_theta.legend(loc="upper right")

    # Number of distinct fixed points Diagram
    P_zetas = [p.Pzeta0NU for p in particles]
    ax_ndfp.set_ylabel("Number of Fixed Points")
    ax_ndfp.scatter(P_zetas, num_of_XP, s=2, color="#E65100", label="X points")
    ax_ndfp.scatter(P_zetas, num_of_OP, s=2, label="O points")
    ax_ndfp.legend()

    plt.ion()
    plt.show(block=True)
