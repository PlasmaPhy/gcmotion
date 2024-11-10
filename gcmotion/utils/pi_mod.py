"""
[0,2π] or [-π,π] mod
--------------------

Simple function to mod an array between [0,2π] or [-π,π].
"""

import numpy as np

from gcmotion.configuration.plot_parameters import energy_contour as config


def pi_mod(x, lim: list = [-np.pi, np.pi]) -> tuple[np.ndarray, list]:
    """Mods an array between [0,2π] or [-π,π].

    .. note::
        array limits must be expressed as multiples of
        ``np.pi``.

    Parameters
    ----------
    x : np.ndarray
        The array to be modded.
    lim : list
        The mod limits.

    Returns
    -------
    2-tuple containing an np.ndarray and a list
        The modded array and list containing the limits.

    """
    if (lim != [0, 2 * np.pi]) and (lim != [-np.pi, np.pi]):
        print("x_lim must be either [0,2*np.pi] or [-np.pi,np.pi].")
        print("Defaulting to [-π, π].")
        return pi_mod(x, lim=[-np.pi, np.pi])

    if lim == [0, 2 * np.pi]:
        x_plot = np.mod(x, 2 * np.pi)
    elif lim == [-np.pi, np.pi]:
        x_plot = np.mod(x, 2 * np.pi)
        x_plot = x_plot - 2 * np.pi * (x_plot > np.pi)
    return x_plot, lim


def yspan(x: np.array) -> tuple(tuple, tuple):

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
