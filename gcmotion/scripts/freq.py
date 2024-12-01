import numpy as np
import matplotlib.pyplot as plt

from gcmotion.entities.profile import Profile
from gcmotion.plot._base._base_profile_contour import _base_profile_contour
from gcmotion.plot._base._base_contour_colorbar import _base_contour_colorbar

from gcmotion.configuration.scripts_configuration import (
    FrequencyConfig as config,
)


def frequency(profile: Profile, **args):
    r"""Finds the :math:`\omega(E)` relation graphically"""

    # Unpack arguements and pop the ones that must be fixed
    args.pop("thetalim", None)
    args["psilim"] = args.get("psilim", config.psilim)
    args["levels"] = args.get("levels", config.levels)
    args["flux_units"] = args.get("flux_units", config.flux_units)
    args["E_units"] = args.get("E_units", config.E_units)
    args["potential"] = args.get("potential", config.potential)
    args["wall"] = args.get("wall", config.wall)
    args["grid_density"] = args.get("grid_density", config.grid_density)
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
    # Create the following contours:
    # 1: Between [-π, π], and keep all the closed lines
    # 2: Between [0, π], and keep all the closed lines
    # 3: Between [-2π, 2π], and keep all the open lines
    # That way we capture every orbit line without having to process them
    # further
    C = _base_profile_contour(
        profile=profile,
        ax=ax_dict["contour"],
        thetalim=[-2 * np.pi, 2 * np.pi],
        **args,
    )
    composite_paths = C.get_paths()

    # NOTE: It seems that the number 79 appears inside path.codes once for
    # every different path.
    paths = split_composite_paths(composite_paths)

    # NOTE: Every contour path that touches the wall (e.g passing) contains
    # both -np.pi and np.pi in its vertices. Classify them accordingly:
    trapped_paths = []
    passing_paths = []
    for path in paths:
        if ispassing(path):
            passing_paths.append(path)
        else:
            trapped_paths.append(path)

    # NOTE: For the passing, make sure that they indeed span (-2π,2π)
    iterpaths = passing_paths
    passing_paths = []
    for path in iterpaths:
        if abs(np.max(path[0]) - np.min(path[0])) > 12:  # 4π = 12.5...
            passing_paths.append(path)

    # NOTE: For all the paths, we must remove the edges and add the two
    # points (-2*np.pi, 0) and (2*np.pi, 0).
    psilim = profile.Q(args["psilim"], "psi_wall").to(args["flux_units"]).m
    cropped_trapped_paths = []
    cropped_passing_paths = []
    for path in trapped_paths:
        cropped_trapped_paths.append(remove_edges(path, psilim))
    for path in passing_paths:
        cropped_passing_paths.append(remove_edges(path, psilim))
    trapped_paths = cropped_trapped_paths
    passing_paths = cropped_passing_paths

    # For the trapped paths, keep the ones that are closed
    iterpaths = trapped_paths
    trapped_paths = []
    for path in iterpaths:
        if path.shape[1] > 4:
            if abs(np.max(np.diff(path, axis=1))) < 0.3:
                trapped_paths.append(path)

    # For all the paths, discard those with points to close to the upper
    # and lower limit
    iterpaths = trapped_paths
    trapped_paths = []
    for path in iterpaths:
        psimax = np.max(path[1])
        if (
            psimax < 0.95 * psilim[1]  # if it touches the wall
            and psimax > 0.03 * psilim[1]  # Too low
        ):
            trapped_paths.append(path)
    iterpaths = passing_paths
    passing_paths = []
    for path in iterpaths:
        psimax = np.max(path[1])
        if psimax < 0.95 * psilim[1]:  # if it touches the wall
            passing_paths.append(path)

    # for all passing paths, add (-2π,0) and (2π,0)
    iterpaths = passing_paths
    passing_paths = []
    for path in iterpaths:
        passing_paths.append(
            np.append(path, [[-2 * np.pi, 2 * np.pi], [0, 0]], axis=1)
        )

    # plot found paths:
    for path in passing_paths:
        ax_dict["passing"].scatter(path[0], path[1], s=1)
        ax_dict["passing"].set_title("Passing")
        ax_dict["passing"].set_xlim([-2 * np.pi, 2 * np.pi])
        ax_dict["passing"].set_ylim(psilim)

    for path in trapped_paths:
        ax_dict["trapped"].scatter(path[0], path[1], s=1)
        ax_dict["trapped"].set_title("Trapped")

    # # Plot them all nicely together
    # _ = _base_profile_contour(
    #     profile=profile,
    #     ax=ax_dict["all"],
    #     thetalim=[-2 * np.pi, 2 * np.pi],
    #     **args,
    # )
    # for path in trapped_paths:
    #     ax_dict["all"].scatter(path[0], path[1], s=1)
    #
    # for path in passing_paths:
    #     ax_dict["all"].scatter(path[0], path[1], s=1)

    # Calculate areas
    areas = []
    energies = []
    for path in trapped_paths + passing_paths:
        x = path[0]
        y = path[1]
        area = 0.5 * np.array(
            np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))
        )
        psi = profile.Q(y[4], args["flux_units"])
        theta = x[4]
        energy = profile.findEnergy(psi, theta, args["E_units"])
        areas.append(area)
        energies.append(energy)

    # Plot relation
    ax_dict["freq"].scatter(energies, areas, s=4)
    ax_dict["freq"].margins(0.01)
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
    return trapped_paths, passing_paths


def split_composite_paths(paths):
    r"""Returns a list of splitted paths"""

    # Iterate over all paths, find the indeces where path.codes == 79 and split
    # accordingly
    path_list = []
    for path in paths:
        idx = np.argwhere(path.codes == 79).flatten()
        extracted_paths = np.split(path.vertices, idx)
        # Trim the paths since the first and last points are random
        # Also removed the last extracted path since it is empty
        # Also transpose them, easier to work with
        path_list += [trimmed[1:-1].T for trimmed in extracted_paths[:-1]]

    return path_list


def ispassing(path):

    thetas = path[0]
    if -2 * np.pi in thetas or 2 * np.pi in thetas:
        return True
    else:
        return False


def remove_edges(path, psilim):

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

    return path
