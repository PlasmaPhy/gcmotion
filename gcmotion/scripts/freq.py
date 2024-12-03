import numpy as np
import matplotlib.pyplot as plt
from math import isclose

from gcmotion.entities.profile import Profile
from gcmotion.plot._base._base_profile_contour import _base_profile_contour

from gcmotion.configuration.scripts_configuration import FrequencyConfig


def frequency(profile: Profile, **args):
    r"""Finds the :math:`\omega(E)` relation graphically"""

    # Unpack parameters
    args.pop("thetalim", None)
    config = FrequencyConfig()
    for key, value in args.items():
        setattr(config, key, value)

    # Create figure and axes
    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
    }
    fig = plt.figure(**fig_kw)
    ax_dict = fig.subplot_mosaic(
        [
            ["contour", "freq"],
            ["passing", "trapped"],
        ]
    )

    # Create the contour from which we extract the paths
    C = _base_profile_contour(
        profile=profile,
        ax=ax_dict["contour"],
        thetalim=[-2 * np.pi, 2 * np.pi],
        **args,
    )
    composite_paths = C.get_paths()

    # It seems that the number 79 appears inside path.codes once for every
    # different path.
    # This also transposes the path vertices
    paths = split_composite_paths(composite_paths)

    # Remove ALL points that lie on the bounding box
    psilim = profile.Q(config.psilim, "psi_wall").to(config.flux_units).m
    paths = remove_edges(paths, psilim)

    # Sort paths' points with respect to their theta coord
    paths = sort_paths(paths)

    for path in paths:
        # ax_dict["passing"].scatter(path[0], path[1], s=1)
        ax_dict["passing"].plot(path[0], path[1])

    plt.show()
    return paths

    # Every contour path that touches the wall (e.g passing) contains both
    # -np.pi and np.pi in its vertices. Classify them accordingly:.
    passing, trapped = classify_paths(paths)

    # # NOTE: For the passing, make sure that they indeed span (-2π,2π)
    # iterpaths = passing
    # passing = []
    # for path in iterpaths:
    #     if abs(np.max(path[0]) - np.min(path[0])) > 12:  # 4π = 12.5...
    #         passing.append(path)
    #
    # # NOTE: For all the paths, we must remove the edges and add the two
    # # points (-2*np.pi, 0) and (2*np.pi, 0).
    # psilim = profile.Q(args["psilim"], "psi_wall").to(args["flux_units"]).m
    # cropped_trapped = []
    # cropped_passing = []
    # for path in trapped:
    #     cropped_trapped.append(remove_edges(path, psilim))
    # for path in passing:
    #     cropped_passing.append(remove_edges(path, psilim))
    # trapped = cropped_trapped
    # passing = cropped_passing
    #
    # # For the trapped paths, keep the ones that are closed
    # iterpaths = trapped
    # trapped = []
    # for path in iterpaths:
    #     if path.shape[1] > 4:
    #         if abs(np.max(np.diff(path, axis=1))) < 0.3:
    #             trapped.append(path)
    #
    # # For all the paths, discard those with points to close to the upper
    # # and lower limit
    # iterpaths = trapped
    # trapped = []
    # for path in iterpaths:
    #     psimax = np.max(path[1])
    #     if (
    #         psimax < 0.95 * psilim[1]  # if it touches the wall
    #         and psimax > 0.03 * psilim[1]  # Too low
    #     ):
    #         trapped.append(path)
    # iterpaths = passing
    # passing = []
    # for path in iterpaths:
    #     psimax = np.max(path[1])
    #     if psimax < 0.95 * psilim[1]:  # if it touches the wall
    #         passing.append(path)
    #
    # # for all passing paths, add (-2π,0) and (2π,0)
    # iterpaths = passing
    # passing = []
    # for path in iterpaths:
    #     passing.append(
    #         np.append(path, [[-2 * np.pi, 2 * np.pi], [0, 0]], axis=1)
    #     )
    #
    # plot found paths:
    for path in passing:
        # ax_dict["passing"].scatter(path[0], path[1], s=1)
        ax_dict["passing"].plot(path[0], path[1])
        ax_dict["passing"].set_title("Passing")
        ax_dict["passing"].set_xlim([-2 * np.pi, 2 * np.pi])
        ax_dict["passing"].set_ylim(psilim)

    for path in trapped:
        # ax_dict["trapped"].scatter(path[0], path[1], s=1)
        ax_dict["trapped"].plot(path[0], path[1])
        ax_dict["trapped"].set_title("Trapped")

    # # Plot them all nicely together
    # _ = _base_profile_contour(
    #     profile=profile,
    #     ax=ax_dict["all"],
    #     thetalim=[-2 * np.pi, 2 * np.pi],
    #     **args,
    # )
    # for path in trapped:
    #     ax_dict["all"].scatter(path[0], path[1], s=1)
    #
    # for path in passing:
    #     ax_dict["all"].scatter(path[0], path[1], s=1)

    # # Calculate areas
    # areas = []
    # energies = []
    # for path in trapped + passing:
    #     x = path[0]
    #     y = path[1]
    #     area = 0.5 * np.array(
    #         np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))
    #     )
    #     psi = profile.Q(y[4], config.flux_units)
    #     theta = x[4]
    #     energy = profile.findEnergy(psi, theta, config.E_units)
    #     area = abs(area)
    #     areas.append(area)
    #     energies.append(energy)
    #
    # Plot relation
    # ax_dict["freq"].scatter(energies, areas, s=4)
    # ax_dict["freq"].semilogx()
    # ax_dict["freq"].margins(0.01)
    # ax_dict["freq"].semilogx()

    # Plot low-level contour
    ax_dict["contour"].clear()
    args["levels"] = 30
    C = _base_profile_contour(
        profile=profile,
        ax=ax_dict["contour"],
        thetalim=[-np.pi, np.pi],
        **args,
    )

    plt.show()
    return passing, trapped


def split_composite_paths(paths):
    r"""Returns a list of splitted paths"""

    # Iterate over all paths, find the indeces where path.codes == 79 and split
    # accordingly
    new_paths = []
    for path in paths:
        idx = np.argwhere(path.codes == 79).flatten()
        extracted = np.split(path.vertices, idx)
        # Trim the paths since the first and last points are random
        # Also removed the last extracted path since it is empty
        # Also transpose them, easier to work with
        new_paths += [trimmed[1:-2].T for trimmed in extracted[:-1]]

    return new_paths


def remove_edges(paths, psilim):

    for n, path in enumerate(paths):
        tol = abs(psilim[1] - psilim[0]) / 1000
        cond = (
            (path[0] > -2 * np.pi)
            & (path[0] < 2 * np.pi)
            & (abs(path[1] - psilim[0]) > tol)
            & (abs(path[1] - psilim[1]) > tol)
        )
        thetas, psis = np.array(np.where(cond, path, np.nan))
        path = np.array([thetas, psis])
        path = path[:, ~np.isnan(path).any(axis=0)]
        paths[n] = path

    return paths


def sort_paths(paths):
    for n, path in enumerate(paths):
        paths[n] = path[:, path[0, :].argsort()]
        del path
    return paths


def classify_paths(paths):
    passing = []
    trapped = []
    for path in paths:
        if path[0, -1] - path[0, 0] > 12:
            passing.append(path)
        else:
            trapped.append(path)
    return passing, trapped
