import numpy as np
import matplotlib.pyplot as plt

from gcmotion.plotters.drift import drift


def collection_drift(collection, angle, lim, **params):

    params["canvas"] = params.get("canvas", None)

    if params["canvas"] is None:
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        params["canvas"] = (fig, ax)
    else:  # Use external canvas
        fig, ax = params["canvas"]

    params["different_colors"] = params.get("different_colors", False)
    params["plot_initial"] = params.get("plot_initial", True)

    for p in collection.particles:
        drift(p, angle, lim, **params)

    fig.set_tight_layout(True)
    plt.ion()
    plt.show(block=True)
