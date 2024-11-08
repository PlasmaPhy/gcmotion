import matplotlib.pyplot as plt
import numpy as np

from gcmotion.scripts.bifurcation import bifurcation
from gcmotion.classes.collection import Collection


def bifurcation_plot(
    collection: Collection,
    theta_lim: list = [-np.pi, np.pi],
    P_theta_lim: list = [0.01, 1.3],
    iterations: int = 250,
    info: bool = False,
):

    thetas_fixed, P_thetas_fixed, num_of_fp = bifurcation(
        collection=collection,
        theta_lim=theta_lim,
        P_theta_lim=P_theta_lim,
        iterations=iterations,
        info=info,
    )

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
    ax_P_theta.set_ylabel(r"$P_{\theta_s}/P_{\theta_w}$")
    ax_P_theta.scatter(P_zeta_plot, P_theta_plot, s=2)
    ax_P_theta.axhline(y=1, color="black", linestyle="-", linewidth=2)

    # Number of distinct fixed points Diagram
    # ax_ndfp.set_title("Number of Fixed Points Bifurcation Diagram")
    P_zetas = [p.Pzeta0 for p in particles]
    ax_ndfp.set_ylabel("Fixed Points")
    ax_ndfp.scatter(P_zetas, num_of_fp, s=2)

    plt.ion()
    plt.show(block=True)
