"""
Simple function that determines which elements [,] of a list of lists of len 2 [[,],[,],[,]...]
can be considered distinct from one another. In this project's context, it is used to 
determine which points [x,y] can be considered distinct.

    Parameters
    ----------
    x : list
        List that contains sublists that are to be examined for uniquness.
    tol : list
        If two sublists have both elements that are less than tol (tolerance) apart, they
        are considered idenical.

    Returns
    -------
    List that contains only the distinct sublists/points.


"""

import numpy as np


def distinctify(points: np.ndarray, tol: float):

    def are_considered_equal(sublist1, sublist2, tol=tol):
        return abs(sublist1[0] - sublist2[0]) <= tol and abs(sublist1[1] - sublist2[1]) <= tol

    distinct_points = np.empty((points.shape[0], points.shape[1]))
    distinct_points[:] = np.nan

    idx = 0

    for point in points:

        is_unique = True
        for distinct in distinct_points:
            if are_considered_equal(distinct, point):
                is_unique = False
                break

        if is_unique:
            distinct_points[idx] = point
            idx += 1

    distinct_points = distinct_points[~np.isnan(distinct_points).any(axis=1)]

    return distinct_points
