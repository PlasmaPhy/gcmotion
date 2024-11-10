import numpy as np

from gcmotion.configuration.plot_parameters import energy_contour as config


def yspan(x: np.array):

    zoomout = config["auto_yaxis_zoom"]
    hardylim = config["hardylim"]

    xmin = x.min()
    xmax = x.max()
    minpos = x.argmin()
    maxpos = x.argmax()

    diff = xmax - xmin
    mid = (xmin + xmax) / 2

    lower = max(0, mid - zoomout * diff)
    higher = min(mid + zoomout * diff, hardylim)

    yspan = (lower, higher)
    minmaxpos = (minpos, maxpos)

    return yspan, minmaxpos
