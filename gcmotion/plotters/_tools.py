import numpy as np

from gcmotion.configuration.plot_parameters import energy_contour as config

from matplotlib.axes._base import _AxesBase

type format_coord = _AxesBase.format_coord


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


def yspan(x: np.ndarray, psi_wall: float) -> tuple[float, float]:
    """Calculates the plot span of a bounded array to be plotted.

    The "zoom out" and "hardylim" factors of the plot span can
    be tweaked in the plot_parameters.py under the name
    "auto_yaxis_zoom".

    Parameters
    ----------
    x : np.ndarray
        The array to be plotted.
    psi_wall : float
        The tokamak's psi_wall. Needed to
        calculate the correct hard-y limit.

    Returns
    -------
    2tuple
        2tuple containing the plot limits.

    """

    zoomout = config["auto_yaxis_zoom"]
    hardylim = config["hardylim"] * psi_wall

    xmin = x.min()
    xmax = x.max()

    diff = xmax - xmin
    mid = (xmin + xmax) / 2

    lower = max(0, mid - zoomout * diff)
    higher = min(mid + zoomout * diff, hardylim)

    xspan = (float(lower), float(higher))

    return xspan


def grid_cursor_format_coord(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    z_label: str = "",
) -> format_coord:
    r"""
    Creates a data cursor on contour plots.

    Parameters
    ----------
    x : np.ndarray
        2D numpy array containing the x grid values.
    y : np.ndarray
        2D numpy array containing the x grid values.
    z : np.ndarray
        2D numpy array containing the contour values.
    z_label : str, optional
        The label next to the z values. Defaults to "".
    """

    Xflat, Yflat, Zflat = (
        x.flatten(),
        y.flatten(),
        z.flatten(),
    )

    def fmt(x, y):
        # get closest point with known data
        dist = np.linalg.norm(np.vstack([Xflat - x, Yflat - y]), axis=0)
        idx = np.argmin(dist)
        z = Zflat[idx]
        return f"theta={x:.4g}  Ptheta={y:.5f}  E={z:.4g}[{z_label}]"

    return fmt
