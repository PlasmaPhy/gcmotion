import numpy as np
import matplotlib.pyplot as plt
from math import isclose
from typing import Callable

from gcmotion.entities.profile import Profile
from gcmotion.plot._base._base_profile_contour import _base_profile_contourE
from gcmotion.plot._base._base_contour_colorbar import _base_contour_colorbar

from gcmotion.configuration.scripts_configuration import FrequencyConfig


def frequencies(profile: Profile, **args):

    # Setup configuration
    args.pop("thetalim", None)
    args.pop("mode", None)
    config = FrequencyConfig()
    for key, value in args.items():
        setattr(config, key, value)

    # Create a phony figure to extract the segments and delete it
    # NOTE: We must set up the energy levels manual in order to use them in
    # both ωθ and ωζ, so we need to find the minimum and maximum energies.
    phony_fig, phony_ax = plt.subplots()
    Cphony = _base_profile_contourE(
        profile=profile,
        ax=phony_ax,
        thetalim=[-2 * np.pi, 2 * np.pi],
        mode="lines",
        **args,
    )
    del phony_fig, phony_ax
    plt.close()

    contours = Cphony.allsegs
    _psilim = [Cphony._mins[1], Cphony._maxs[1]]
    Emin = Cphony.zmin
    Emax = Cphony.zmax
    Espan = np.linspace(Emin, Emax, config.levels)

    result = freq_from_contours(
        contours,
        _psilim,
        Espan,
        profile.findEnergy,
        config.flux_units,
        config.E_units,
        profile.Q,
    )
    passing_frequencies = result["passing_frequencies"]
    passing_energies = result["passing_energies"]
    trapped_frequencies = result["trapped_frequencies"]
    trapped_energies = result["trapped_energies"]

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

    # "freq" ax exists in all cases, so start with that
    ax_dict["freq"].scatter(
        passing_energies, passing_frequencies, c="r", marker=".", s=20
    )
    ax_dict["freq"].scatter(
        trapped_energies, trapped_frequencies, c="b", marker=".", s=20
    )
    ax_dict["freq"].loglog()

    # Plot low-level contour
    if "contour" in ax_dict.keys():
        ax_dict["contour"].clear()
        args["levels"] = 30
        Cnew = _base_profile_contourE(
            profile=profile,
            ax=ax_dict["contour"],
            thetalim=[-np.pi, np.pi],
            mode="filled",
            **args,
        )

        cbar = fig.colorbar(Cnew, cax=None, ax=ax_dict["contour"])
        _base_contour_colorbar(ax=cbar.ax, contour=Cnew, numticks=10)
        cbar.ax.set_title(label=f"Energy [{config.E_units}]", size=10)

    if config.mosaic == "debug":
        for path in result["passing_paths"]:
            ax_dict["passing"].plot(path[0], path[1])
            ax_dict["passing"].set_xlim([-2 * np.pi, 2 * np.pi])
            ax_dict["passing"].set_ylim(_psilim)
        for path in result["trapped_paths"]:
            ax_dict["trapped"].plot(path[0], path[1])

    plt.show()


def freq_from_contours(
    segments, psispan, Espan, findEnergy, flux_units, E_units, Q
):
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
    passing_areas = Q(passing_areas, f"{E_units}*seconds")
    trapped_areas = Q(trapped_areas, f"{E_units}*seconds")

    # For every path, calculate its corresponding energy using the first (θ,ψ)
    # pair. Its easier and safer than reusing the Espan, since we cant be sure
    # that the order is the same, and the performance hit is negligible.
    passing_energies = path_energies(
        passing, findEnergy, flux_units, E_units, Q
    )
    trapped_energies = path_energies(
        trapped, findEnergy, flux_units, E_units, Q
    )

    # Calculate energies from ω = dH/dJ, where J is the area we calculated and
    # dH is the approximate derivative of the hamiltonian we calculated. The
    # approximation is better the more levels we take.
    try:
        passing_frequencies = abs(np.gradient(passing_energies, passing_areas))
    except IndexError:
        passing_frequencies = []
    try:
        trapped_frequencies = abs(np.gradient(trapped_energies, trapped_areas))
    except IndexError:
        trapped_frequencies = []

    return {
        "passing_paths": passing,
        "passing_frequencies": passing_frequencies,
        "passing_energies": passing_energies,
        "trapped_paths": trapped,
        "trapped_frequencies": trapped_frequencies,
        "trapped_energies": trapped_energies,
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
