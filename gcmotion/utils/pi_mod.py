"""
[0,2π] or [-π,π] mod
--------------------

Simple function to mod an array between [0,2π] or [-π,π].
"""

import numpy as np


def pi_mod(x, x_lim=[-np.pi, np.pi]):
    """Mods an array between [0,2π] or [-π,π].

    Parameters
    ----------
    x : np.ndarray
        The array to be modded.
    x_lim : list
        The mod limits.

    Returns
    -------
    np.ndarray
        The modded array.

    """
    if (x_lim != [0, 2 * np.pi]) and (x_lim != [-np.pi, np.pi]):
        print("x_lim must be either [0,2*np.pi] or [-np.pi,np.pi].")
        print("Defaulting to [-π, π].")
        return pi_mod(x, x_lim=[-np.pi, np.pi])

    if x_lim == [0, 2 * np.pi]:
        x_plot = np.mod(x, 2 * np.pi)
    elif x_lim == [-np.pi, np.pi]:
        x_plot = np.mod(x, 2 * np.pi)
        x_plot = x_plot - 2 * np.pi * (x_plot > np.pi)
    return x_plot, x_lim
