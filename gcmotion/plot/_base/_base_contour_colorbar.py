import numpy as np
from matplotlib.ticker import MaxNLocator

from matplotlib.axes import Axes
from matplotlib.contour import QuadContourSet

from gcmotion.plot._base._config import _ColorbarConfig


def _base_contour_colorbar(ax: Axes, contour: QuadContourSet, **kwargs):
    r"""Customizes an *existing* colorbar.

    The colorbar must have been already created on a figure, and its reference
    axis is to be passed to this function.

    Notes
    -----
    For a full list of all available optional parameters, see the dataclass
    _ProfileContour at gcmotion/plot/_base/_config. The defaults values are set
    there, and are overwritten if passed as arguements.
    """

    # Unpack parameters
    config = _ColorbarConfig()
    for key, value in kwargs.items():
        setattr(config, key, value)

    # Get contour span
    # Dont use the absolute min and max because it creates gaps in the edges of
    # the colorbar
    levels = np.sort(contour.levels)
    span = [levels[0], levels[-1]]

    # Set tick values
    locator = MaxNLocator(nbins=config.numticks - 1)
    ticks = locator.tick_values(span[0], span[-1])
    # Occasionally, ticks fall outside the energy span and create gaps in the
    # colobar, so we need to fix it
    ticks = np.clip(ticks, levels[0], levels[-1])
    ticks = np.unique(ticks)

    # Draw ticks on colorbar's ax
    ax.set_yticks(
        ticks=ticks,
        labels=[f"{x:.4g}" for x in ticks],
    )

    # Set the passed label
    ax.set_title(label=config.label, size=config.labelsize)
