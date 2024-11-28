import matplotlib.pyplot as plt

from gcmotion.plot.base.base_profile_contour import base_profile_contour
from gcmotion.configuration.plot_parameters import (
    ProfileContourConfig as config
)


def profile_contour():

    # Create figure
    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
    }
    fig = plt.figure(**fig_kw)
    contour = fig.suplots()
