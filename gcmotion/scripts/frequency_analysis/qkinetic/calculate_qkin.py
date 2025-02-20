from gcmotion.configuration.scripts_configuration import CalculateQkinConfig
from gcmotion.scripts.frequency_analysis.contours.contour_generators import (
    main_contour,
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
    profile.psilim = psilim
    profile.valid_orbits = valid_orbits

    debug_plot_valid_orbits(profile)

    return valid_orbits
