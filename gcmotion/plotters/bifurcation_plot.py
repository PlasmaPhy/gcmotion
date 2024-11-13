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
    theta_density=5,
    P_theta_density=5,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: list = [0.01, 1.3],
    info: bool = False,
):
    r"""Draws the bifurcation diagrams for the :math:`theta`'s  fixed,
    the :math:`P_{theta}`'s fixed and the number of fixed points found for
    each :math:`P_{\zeta0}`.

    :meta public:

    Parameters
    ----------
    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles
    theta_density : int, optional
        Integer dictating the number of initial conditions with regard to the
        :math:`\theta` variable that will be passed into :py:func:`bifurcation`
        and then :py:func:`fixed_points`.
    P_theta_density : int, optional
        Integer dictating the number of initial conditions with regard to the
        :math:`P_{\theta}` variable that will be passed into :py:func:`bifurcation`
        and then :py:func:`fixed_points`.
    theta_lim : list, optional
        Provides the limits for the solution search area for fixed points
        with regards to the :math:`\theta` variable. It will be passed into
        :py:func:`bifurcation`
    psi_lim : list, optional
        Provides the limits (divided by psi_wall) for the solution search area with regards
        to the :math:`P_{\theta}` variable. It will be passed into the "bounds" argument of
        :py:func:`bifurcation`.
    info : bool, optional
        Boolean that dictates weather the :math:`P_{\zeta0}` of the particle whose
        fixed points have just been calculated, will be printed alongside the fixed points
        found.
    """

    start = time()
    thetas_fixed, P_thetas_fixed, num_of_fp = bifurcation(
        collection=collection,
        theta_density=theta_density,
        P_theta_density=P_theta_density,
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        info=info,
    )

    print(f"BIFURCATION RUN IN {(time() - start)/60:.1f} mins")

    particles = collection.particles
    p1 = particles[0]
    psi_wall = p1.psi_wall

    fig, ax = plt.subplots(3, 1, figsize=(6, 4), sharex=True)
    plt.xlabel(r"$P_{\zeta}$")

    ax_theta = ax[0]
    ax_P_theta = ax[1]
    ax_ndfp = ax[2]

    P_zeta_plot = []
    theta_plot = []
    P_theta_plot = []

    particles = collection.particles

    # Theta Bifurcation
    for i, p in enumerate(particles):
        P_zeta = p.Pzeta0
        y_list = thetas_fixed[i]
        P_zeta_plot.extend([P_zeta] * y_list.shape[0])
        theta_plot.extend(y_list)

    # ax_theta.set_title("theta Bifurcation Diagram")
    ax_theta.set_ylabel(r"$\theta_s$")
    ax_theta.set_yticks([-np.pi, 0, np.pi])
    ax_theta.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax_theta.set_ylim([-np.pi - 0.5, np.pi + 0.5])
    ax_theta.scatter(P_zeta_plot, theta_plot, s=2)

    P_zeta_plot = []

    # P_theta Bifurcation
    for i, p in enumerate(particles):
        P_zeta = p.Pzeta0
        y_list = P_thetas_fixed[i]
        P_zeta_plot.extend([P_zeta] * y_list.shape[0])
        P_theta_plot.extend(y_list)

    P_theta_plot = [P_theta / psi_wall for P_theta in P_theta_plot]

    # ax_P_theta.set_title(r"$P_{theta}$ Bifurcation Diagram")
    ax_P_theta.set_ylabel(r"$P_{\theta_s}$")
    ax_P_theta.scatter(P_zeta_plot, P_theta_plot, s=2)
    ax_P_theta.axhline(y=1, color="black", linestyle="--", linewidth=1, alpha=0.5)

    # Number of distinct fixed points Diagram
    # ax_ndfp.set_title("Number of Fixed Points Bifurcation Diagram")
    P_zetas = [p.Pzeta0 for p in particles]
    ax_ndfp.set_ylabel("Fixed Points")
    ax_ndfp.scatter(P_zetas, num_of_fp, s=2)

    plt.ion()
    plt.show(block=True)
