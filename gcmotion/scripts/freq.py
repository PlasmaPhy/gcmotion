import numpy as np
import matplotlib.pyplot as plt
from math import isclose

from gcmotion.entities.profile import Profile
from gcmotion.plot._base._base_profile_contour import _base_profile_contour
from gcmotion.plot._base._base_contour_colorbar import _base_contour_colorbar

from gcmotion.configuration.scripts_configuration import FrequencyConfig


def frequency(profile: Profile, **args):

    config = FrequencyConfig()
    for key, value in args.items():
        setattr(config, key, value)

    # Create a phony figure to extract the segments and delete it
    phony_fig, phony_ax = plt.subplots()
    C = _base_profile_contour(
        profile=profile,
        ax=phony_ax,
        thetalim=[-2 * np.pi, 2 * np.pi],
        mode="lines",
        **args,
    )
    result = _analyze_segments(C.allsegs, **args)
    del phony_fig, phony_ax

    passing = result["passing"]
    trapped = result["trapped"]

    # Create figure and axes
    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
    }
    fig = plt.figure(**fig_kw)

    match config.mosaic:
        case "simple":
            mosaic = [["contour", "freq"]]
        case "freq":
            mosaic = [["freq"]]
        case _:
            mosaic = [
                ["contour", "freq"],
                ["passing", "trapped"],
            ]
    ax_dict = fig.subplot_mosaic(mosaic)


def _analyze_segments(segments, **args):
    r"""Finds the :math:`\omega(E)` relation graphically"""

    # Unpack parameters
    args.pop("thetalim", None)
    args.pop("mode", None)
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
        mode="lines",
        **args,
    )
    segments = C.allsegs

    # Unpack segments and store as paths (2xN)
    paths = flatten_segments(segments)

    # Remove empty paths
    paths = remove_empty(paths)

    # Sort paths' points with respect to their theta coord
    # paths = sort_paths(paths)

    # Every contour path that touches the wall (e.g passing) contains both
    # -2np.pi and 2np.pi in its vertices. Classify them accordingly:.
    passing, trapped = classify_paths(paths)

    # Paths that get cut off by the bounding box are classified as trapped, so
    # remove them.
    _psilim = profile.Q(config.psilim, "psi_wall").to(config.flux_units).m
    trapped = remove_cutoffs(trapped, _psilim)

    passingE, trappedE = shoelace(passing, trapped)

    ############
    # Plotting #
    ############

    # Passing/Trapped contours
    for path in passing:
        ax_dict["passing"].set_xlim([-2 * np.pi, 2 * np.pi])
        ax_dict["passing"].set_ylim(_psilim)
        ax_dict["passing"].plot(path[0], path[-1])
        if config.st_end_points:
            ax_dict["passing"].scatter(
                path[0, 0], path[1, 0], marker="+", s=20
            )
            ax_dict["passing"].scatter(
                path[0, -1], path[1, -1], marker="x", s=20
            )

    for path in trapped:
        ax_dict["trapped"].plot(path[0], path[-1])
        if config.st_end_points:
            ax_dict["trapped"].scatter(
                path[0, 0], path[1, 0], marker="+", s=20
            )
            ax_dict["trapped"].scatter(
                path[0, -1], path[1, -1], marker="x", s=20
            )

    # Areas
    for n, path in enumerate(passing):
        psi = profile.Q(path[1, 0], config.flux_units)
        theta = path[0, 0]
        Energy = profile.findEnergy(psi, theta, config.E_units)
        ax_dict["freq"].scatter(Energy, passingE[n], marker=".", s=20, c="red")

    for n, path in enumerate(trapped):
        psi = profile.Q(path[1, 0], config.flux_units)
        theta = path[0, 0]
        Energy = profile.findEnergy(psi, theta, config.E_units)
        ax_dict["freq"].scatter(
            Energy, trappedE[n], marker=".", s=20, c="blue"
        )
    ax_dict["freq"].semilogx()
    ax_dict["freq"].margins(0.02)

    ax_dict["passing"].clear()
    # Plot low-level contour
    ax_dict["contour"].clear()
    args["levels"] = 30
    Cnew = _base_profile_contour(
        profile=profile,
        ax=ax_dict["contour"],
        thetalim=[-np.pi, np.pi],
        mode="filled",
        **args,
    )

    cbar = fig.colorbar(Cnew, cax=None, ax=ax_dict["contour"])
    _base_contour_colorbar(ax=cbar.ax, contour=Cnew, numticks=10)
    cbar.ax.set_title(label=f"Energy [{config.E_units}]", size=10)

    plt.show()
    return C


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
        E = 0.5 * np.array(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
        return abs(E)

    # Add bottom points
    passing_full = add_bottom_points(passing.copy())
    passingE = [0.5 * area(path) for path in passing_full]
    trappedE = [area(path) for path in trapped]

    return passingE, trappedE
