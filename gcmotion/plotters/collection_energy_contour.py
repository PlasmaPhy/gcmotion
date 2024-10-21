import numpy as np
import matplotlib.pyplot as plt

from gcmotion.plotters.collection_drift import collection_drift
from gcmotion.plotters.energy_contour import energy_contour


def collection_energy_contour(
    collection,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: str | list = "auto",
    plot_drift: bool = True,
    contour_Phi: bool = True,
    units: str = "keV",
    levels: int = None,
    wall_shade: bool = True,
    **params,
):

    fig2 = plt.figure(figsize=(6, 4))
    ax2 = fig2.add_subplot(111)
    canvas = (fig2, ax2)
    params["canvas"] = canvas

    # plt.figure()

    energy_contour(
        collection.particles[int(len(collection.particles) / 2)],
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        plot_drift=plot_drift,
        contour_Phi=contour_Phi,
        units=units,
        levels=levels,
        wall_shade=wall_shade,
        **params,
    )

    if plot_drift:
        # params["canvas"] = plt.gcf(), plt.gca()
        collection_drift(collection, "theta", theta_lim, **params)

    fig2.set_tight_layout(True)
    plt.ion()
    plt.show(block=True)
