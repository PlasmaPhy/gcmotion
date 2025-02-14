import numpy as np


def shoelace(array: np.ndarray) -> float:
    r"""Calculates the area of a polygon.

    Parameters
    ----------

    array: np.ndarray
        (N,2) numpy array containing the polygon's points.

    Returns
    -------
    float
        The area of the polygon.

    """
    x, y = array.T[:]
    return float(
        0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    )
