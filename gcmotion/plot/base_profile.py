import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger
from gcmotion.configuration.plot_parameters import (
    BaseProfileContourConfig as config,
)
from gcmotion.entities.profile import Profile


def base_profile_contour(profile: Profile, *args):
    r"""Plots a basic theta-Ptheta-Hamiltonian contour."""

    logger.info("==> Plotting Base Profile Contour...")

    # Unpack parameters
    psilim = args.get("psilim", [0, 1.2])
    Pthetalim = args.get("Pthetalim", None)
    ax = args.get("ax", None)

    # Create figure
    # An ax is created unless an already existing one is passed in args.
    fig = plt.figure(**config.fig_kw)
    ax = ax if ax is not None else fig.subplots()

    # Setup meshgrida
