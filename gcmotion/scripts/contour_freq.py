import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from time import time
from statistics import mean
from dataclasses import asdict
from collections import deque
from matplotlib.contour import QuadContourSet
from matplotlib.patches import Rectangle

from gcmotion.utils.logger_setup import logger
from gcmotion.entities.profile import Profile
from gcmotion.configuration.scripts_configuration import ContourFreqConfig
from gcmotion.scripts.utils.contour_segment import ContourSegment
from gcmotion.plot._base._base_profile_energy_contour import (
    _base_profile_energy_contour,
)


def omegas(profile: Profile, psilim: list, **kwargs):

    C = _contour_set(profile, psilim, **kwargs)

    segs = _filter_paths(C, profile)

    if kwargs.get("show_segments"):
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax = fig.add_subplot()
        fig.set_constrained_layout(True)
        for seg in segs:
            seg.simple_plot(ax=ax, show=False)
        plt.show()

    for seg in tqdm(segs, desc="Calculating omegas"):
        _segment_omega(seg, profile, C.data)

    segs = [seg for seg in segs if hasattr(seg, "omega_theta")]

    return segs


def _segment_omega(
    segment: ContourSegment,
    profile: Profile,
    data: dict,
    **kwargs,
):
    r"""Calculates the ωθ of a single ContourSegment, independently from all
    the others."""
    if not segment.valid:
        return

    rtol = 10e-3
    E = segment.E
    Elower = E * (1 - rtol)
    Eupper = E * (1 + rtol)

    # Create a contour
    C = plt.contour(
        "theta", "ycoord", "Energy", data=data, levels=[Elower, E, Eupper]
    )
    C.psilim = segment.psilim
    segments = _filter_paths(C, profile, pbar=False)

    central_seg = [seg for seg in segments if seg.E == E][0]
    upper_segs = [seg for seg in segments if seg.E > E]
    lower_segs = [seg for seg in segments if seg.E < E]

    upper_distances = [
        segment.bbox_distance((seg.theta_min, seg.psi_min))
        for seg in upper_segs
    ]
    lower_distances = [
        segment.bbox_distance((seg.theta_min, seg.psi_min))
        for seg in lower_segs
    ]

    try:
        correct_upper_seg = upper_segs[np.argmin(upper_distances)]
        correct_lower_seg = lower_segs[np.argmin(lower_distances)]
    except ValueError:
        logger.trace("Skipped segment")
        return

    # Sanity check
    assert central_seg.E == E
    assert correct_lower_seg.E == Elower
    assert correct_upper_seg.E == Eupper

    dE1 = correct_upper_seg.E - central_seg.E
    dE2 = central_seg.E - correct_lower_seg.E
    dE = mean([dE1, dE2])

    central_seg.calculate_area()
    correct_upper_seg.calculate_area()
    correct_lower_seg.calculate_area()

    dJ1 = correct_upper_seg.J - central_seg.J
    dJ2 = central_seg.J - correct_lower_seg.J
    dJ = mean([dJ1, dJ2])

    omega = dE / dJ

    segment.omega_theta = omega


def _contour_set(profile: Profile, psilim: list, **kwargs):

    config = ContourFreqConfig()
    for key, value in kwargs.items():
        setattr(config, key, value)

    # Create phony figure
    fig = plt.figure()
    ax = fig.add_subplot()

    kw = {
        **asdict(config),
        "thetalim": [-2 * np.pi, 2 * np.pi],
        "psilim": psilim,
        "levels": kwargs["levels"],
        "E_units": "NUJoule",
        "ycoord": "Ptheta",
        "flux_units": "NUMagnetic_flux",
        "locator": "log",
        "mode": "lines",
    }
    C = _base_profile_energy_contour(profile, ax=ax, **kw)

    psilim = profile.Q(psilim, "psi_wall").to("NUmagnetic_flux").m
    C.psilim = psilim

    if config.fullshow:
        plt.title("1. Full Contour")
        plt.show()
    else:
        plt.close()

    return C


def _filter_paths(
    C: QuadContourSet, profile: Profile, pbar: bool = True
) -> None:

    contour_segments = deque()

    psilim = C.psilim
    total = sum(len(seg) for seg in C.allsegs)  # Slow but nice progres bar :)
    pbar = tqdm(total=total, desc="Creating ContourSegments", disable=not pbar)

    # Iterate over allsegs and create ContourSegments
    for E, composite_path in zip(C.levels, C.get_paths()):
        # HACK: if this ever changes, just copy the code, its about 10 lines.
        # It should be a public method imho.
        segment_gen = composite_path._iter_connected_components()

        for segment in segment_gen:
            if segment.vertices.shape == (0, 2):
                continue
            newsegment = ContourSegment(
                segment=segment, profile=profile, psilim=psilim, E=E
            )
            contour_segments.append(newsegment)

            pbar.update()
    pbar.close()

    valid_segments = [
        seg for seg in contour_segments if getattr(seg, "valid", False)
    ]
    logger.debug(
        f"Found {len(contour_segments)} total, kept {len(valid_segments)}."
    )

    plt.title("2. Unpacked Paths")

    return valid_segments
