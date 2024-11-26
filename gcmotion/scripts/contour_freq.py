import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.energy_Ptheta import energy_Ptheta


def contour_freq(profile, psi_lim):

    Q = profile.Q

    psi_span = Q(psi_lim, "NUpsi_wall").to("NUMagnetic_flux")
    psi, theta = np.meshgrid(
        np.linspace(psi_span[0].m, psi_span[1].m, 200),
        np.linspace(-np.pi, np.pi, 200),
    )

    W, Ptheta = energy_Ptheta(psi, theta, profile, contour_Phi=True)

    fig, ax = plt.subplots(1, 1)

    C = ax.contourf(theta, Ptheta, W, levels=20)
    plt.ion()
    plt.show(block=True)
    print("BAM")
