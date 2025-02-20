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
        print("No valid orbits found.")
        return

    profile.psilim = psilim
    debug_plot_valid_orbits(profile, valid_orbits)

    # Step 2: Create adjecent Pzeta profiles. Fastest way I could find.
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
        orbit.qkinetic = calculate_orbit_qkin(orbit=orbit, profiles=Profiles)


def calculate_orbit_qkin(orbit: ContourOrbit, profiles: dict) -> float:
    r"""Step 2: For each contour orbit found, calculate qkinetic."""

    config = CalculateQkinConfig()

    # Step 2a: Calculate orbits bounding box
    orbit.calculate_bbox()

    # Step 2b: Create 2 local contours, 1 from each adjacent profile.
    LowerContour = local_contour(profile=profiles["Lower"], bbox=orbit.bbox)
    UpperContour = local_contour(profile=profiles["Upper"], bbox=orbit.bbox)
    profiles["Lower"].psilim = LowerContour["psilim"]
    profiles["Upper"].psilim = UpperContour["psilim"]

    # Step 2c: Try to find contour lines with the same E in each Contour
    lower_orbits = generate_contour_orbits(
        Contour=LowerContour, level=profiles["Main"].ENU.m, config=config
    )
    upper_orbits = generate_contour_orbits(
        Contour=UpperContour, level=profiles["Main"].ENU.m, config=config
    )

    print(f"Lower orbits found: {len(lower_orbits)}")
    print(f"Upper orbits found: {len(upper_orbits)}")

    debug_plot_valid_orbits(profiles["Lower"], lower_orbits)
    debug_plot_valid_orbits(profiles["Upper"], upper_orbits)
