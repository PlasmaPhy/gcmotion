import numpy as np
from copy import deepcopy

from gcmotion.configuration.scripts_configuration import CalculateQkinConfig
from gcmotion.scripts.frequency_analysis.contours.contour_orbits import (
    ContourOrbit,
)
from gcmotion.scripts.frequency_analysis.contours.contour_generators import (
    main_contour,
    local_contour,
)
from gcmotion.scripts.frequency_analysis.contours.lines_processing import (
    generate_contour_orbits,
)

from gcmotion.scripts.frequency_analysis.plots import debug_plot_valid_orbits
from gcmotion.entities.profile import Profile


def calculate_qkin(profile: Profile, psilim: tuple = (0.05, 1), **kwargs):

    config = CalculateQkinConfig()
    for key, value in kwargs.items():
        setattr(config, key, value)

    # Step 1: Generate valid isoenergy orbits. Store them in Profile for
    # conviniency.
    MainContour = main_contour(profile=profile, psilim=psilim)

    valid_orbits = generate_contour_orbits(
        Contour=MainContour, level=profile.ENU.m, config=config
    )

    if len(valid_orbits) == 0:
        # print("No valid orbits found.")
        return [0]

    profile.psilim = psilim
    # debug_plot_valid_orbits(profile, valid_orbits)

    # Step 2: Create adjecent Pzeta profiles. Fastest way I could find.
    profile.psiwall_lim = psilim
    lower_profile, upper_profile = deepcopy(profile), deepcopy(profile)
    lower_profile.PzetaNU = profile.PzetaNU * (1 - config.pzeta_rtol)
    upper_profile.PzetaNU = profile.PzetaNU * (1 + config.pzeta_rtol)
    assert upper_profile is not lower_profile

    profile.valid_orbits = valid_orbits

    Profiles = {
        "Main": profile,
        "Lower": lower_profile,
        "Upper": upper_profile,
    }

    # Step 3: Calculating qkinetics for every contour orbit.
    for orbit in profile.valid_orbits:
        calculate_orbit_qkin(main_orbit=orbit, profiles=Profiles)

    return [orbit for orbit in valid_orbits if hasattr(orbit, "qkinetic")]


def calculate_orbit_qkin(main_orbit: ContourOrbit, profiles: dict) -> float:
    r"""Step 2: For each contour orbit found, calculate qkinetic."""

    config = CalculateQkinConfig()

    # Step 2a: Calculate orbits bounding box and t/p classification
    main_orbit.calculate_bbox()
    main_orbit.classify_as_tp()

    # Step 2b: Create 2 local contours, 1 from each adjacent profile.
    LowerContour = local_contour(profile=profiles["Lower"], orbit=main_orbit)
    UpperContour = local_contour(profile=profiles["Upper"], orbit=main_orbit)
    profiles["Lower"].psilim = LowerContour["psilim"]
    profiles["Upper"].psilim = UpperContour["psilim"]

    # Step 2c: Try to find contour lines with the same E in each Contour
    lower_orbits = generate_contour_orbits(
        Contour=LowerContour, level=profiles["Main"].ENU.m, config=config
    )
    upper_orbits = generate_contour_orbits(
        Contour=UpperContour, level=profiles["Main"].ENU.m, config=config
    )

    # print(f"Lower orbits found: {len(lower_orbits)}")
    # print(f"Upper orbits found: {len(upper_orbits)}")

    # debug_plot_valid_orbits(
    #     profiles["Main"], (*lower_orbits, main_orbit, *upper_orbits)
    # )
    # debug_plot_valid_orbits(profiles["Lower"], lower_orbits)
    # debug_plot_valid_orbits(profiles["Upper"], upper_orbits)

    # Step 2d: Return if orbit is about to disappear
    if len(lower_orbits) == 0 or len(upper_orbits) == 0:
        main_orbit.edge_orbit = True
        return
    else:
        main_orbit.edge_orbit = False

    # Step 2e: Validate orbits, and calculate
    for orbit in lower_orbits:
        orbit.validate(profiles["Lower"].psilim)
        orbit.calculate_bbox()
    for orbit in upper_orbits:
        orbit.validate(profiles["Upper"].psilim)
        orbit.calculate_bbox()

    # Step 2f: pick closest adjacent orbits
    lower_distances = (
        orbit.distance_from(lower_orbit.bbox) for lower_orbit in lower_orbits
    )
    upper_distances = (
        orbit.distance_from(upper_orbit.bbox) for upper_orbit in upper_orbits
    )

    lower_orbit = lower_orbits[np.argmin(lower_distances)]
    upper_orbit = upper_orbits[np.argmin(upper_distances)]

    # Step 3a: Classsify the 2 adjacent orbits as trapped/passing and and base
    # points if needed
    # Step 3b: Converter psis to pthetas
    # Step 3c : Calculate Jtheta
    for orbit, profile in zip(
        [lower_orbit, upper_orbit], [profiles["Lower"], profiles["Upper"]]
    ):
        orbit.classify_as_tp()
        orbit.close_segment()
        orbit.convert_to_ptheta(profile)
        orbit.calculate_Jtheta()

    # Step 3d: Calculate qkinetic
    dJtheta = upper_orbit.Jtheta - lower_orbit.Jtheta
    dPzeta = profiles["Upper"].PzetaNU.m - profiles["Lower"].PzetaNU.m

    qkinetic = -dJtheta / dPzeta
    if abs(qkinetic) > 8:
        return
    else:
        main_orbit.Pzeta = profiles["Main"].PzetaNU.m
        main_orbit.qkinetic = qkinetic
