import numpy as np
from matplotlib.ticker import MaxNLocator

from matplotlib.axes import Axes
from matplotlib.contour import QuadContourSet


def _base_contour_colorbar(ax: Axes, contour: QuadContourSet, **args):
    r"""Customizes an *existing* colorbar.

    The colorbar must have been already created on a figure, and its reference
    axis is to be passed to this function.
    """

    # Unpack parameters
    numticks = args.get("numticks")
    label = args.get("cbarlabel")
    labelsize = args.get("cbarlabelsize")

    # Get contour span
    # Dont use the absolute min and max because it creates gaps in the edges of
    # the colorbar
    levels = np.sort(contour.levels)
    span = [levels[2], levels[-3]]

    # Set tick values
    locator = MaxNLocator(nbins=numticks - 1)
    ticks = locator.tick_values(span[0], span[1])

    # Draw ticks on colorbar's ax
    ax.set_yticks(
        ticks=ticks,
        labels=[f"{x:.4g}" for x in ticks],
    )

    # Set the passed label
    ax.set_title(label=label, size=labelsize)
