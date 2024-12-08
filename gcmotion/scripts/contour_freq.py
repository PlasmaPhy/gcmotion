import numpy as np
from math import isclose
from typing import Callable


def contour_freq(segments, psispan, Espan, findEnergy, flux_units, E_units, Q):
    r"""Finds the :math:`\omega(E)` relation graphically"""

    # Unpack segments and store as paths (2xN)
    paths = flatten_segments(segments)

    # Remove empty paths
    paths = remove_empty(paths)

    # Every contour path that touches the wall (e.g passing) contains both
    # -2np.pi and 2np.pi in its vertices exaclty. Classify them accordingly:
    passing, trapped = classify_paths(paths)

    # Paths that get cut off by the bounding box are classified as trapped, so
    # remove them.
    trapped = remove_cutoffs(trapped, psispan)

    # Calculate the areas with the shoelace algorithm
    passing_areas, trapped_areas = shoelace(passing, trapped)

    return {
        "passing_paths": passing,
        "trapped_paths": trapped,
    }


def flatten_segments(segments):
    paths = []
    for group in segments:
        for segment in group:
            paths.append(np.array(segment).T)
    return paths


def split_composite_paths(paths):
    r"""Returns a list of splitted paths"""

    # Iterate over all paths, find the indeces where path.codes == 79 and split
    # accordingly
    new_paths = []
    for path in paths:
        # print(path)
        idx = np.argwhere(path.codes == 79).flatten()
        # idx += 1
        # print(idx)
        extracted = np.split(path.vertices, idx)
        extracted = [x for x in extracted if x.shape != (1, 2)]
        # print(extracted)
        # Trim the paths since the first and last points are random
        # Also removed the last extracted path since it is empty
        # Also transpose them, easier to work with
        new_paths += [trimmed[:].T for trimmed in extracted[:-1]]

    return new_paths


def remove_empty(paths):
    return [path for path in paths if path.shape[1] > 2]


def sort_paths(paths):
    for n, path in enumerate(paths):
        paths[n] = path[:, path[0, :].argsort()]
        del path
    return paths


def classify_paths(paths):
    passing = []
    trapped = []
    for path in paths:
        if np.max(path[0]) - np.min(path[0]) > 12:
            passing.append(path)
        else:
            trapped.append(path)
    return passing, trapped


def remove_cutoffs(passing, psilim):

    paths = []
    for path in passing:
        thetas = path[0]
        psis = path[1]
        if not (
            isclose(thetas[0], 2 * np.pi)
            or isclose(thetas[0], -2 * np.pi)
            or isclose(thetas[-1], 2 * np.pi)
            or isclose(thetas[-1], -2 * np.pi)
            or isclose(psis[0], psilim[0])
            or isclose(psis[0], psilim[1])
            or isclose(psis[-1], psilim[0])
            or isclose(psis[-1], psilim[1])
        ):
            paths.append(path)
    return paths


def add_bottom_points(passing):
    points = np.array([[-2 * np.pi, 2 * np.pi], [0, 0]])
    for n, path in enumerate(passing):
        passing[n] = np.append(path, points, axis=1)
    return passing


def shoelace(passing, trapped):
    def area(path):
        x = path[0]
        y = path[1]
        area = 0.5 * np.array(
            np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))
        )
        return area

    # Add bottom points
    passing_full = add_bottom_points(passing.copy())
    passingE = [0.5 * area(path) for path in passing_full]
    trappedE = [area(path) for path in trapped]

    return passingE, trappedE


def path_energies(
    paths: np.ndarray, findEnergy: Callable, flux_units: str, E_units: str, Q
):
    if len(paths) == 0:
        return []

    energies = []
    for path in paths:
        psi = Q(path[1, 0], flux_units)
        theta = path[0, 0]
        E = findEnergy(psi, theta, E_units)
        energies.append(E)
    return energies
