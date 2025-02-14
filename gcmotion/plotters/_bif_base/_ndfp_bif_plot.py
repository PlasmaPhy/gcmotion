from collections import deque
from matplotlib.ticker import MaxNLocator
from gcmotion.utils.logger_setup import logger


def _ndfp_bif_plot(
    profiles: list,
    num_of_XP: list | deque,
    num_of_OP: list | deque,
    which_COM: str,
    ax=None,
):
    r"""Base plotting function. Only draws upon a given axis without showing
    any figures.

    Parameters
    ----------
    profiles : list, deque
        List of profile objects.
    num_of_XP : deque, list
        The number of distinct X points found for each COM value.
    num_of_OP : deque, list
        The number of distinct O points found for each COM value.
    which_COM : str
        String that indicates with respect to which constant of motion (COM) :math:`\mu`
        or :math:`P_{\zeta}` the number of distinct fixed points fixed are plotted.
    ax : Axes
        The ax upon which to draw.
    """
    logger.info("\t==> Plotting Base Number of distinct Fixed Points Bifurcation Diagram...")

    selected_COMNU_str = which_COM + "NU"

    # Number of distinct fixed points Diagram
    selected_COMs = [getattr(profile, selected_COMNU_str, "PzetaNU").m for profile in profiles]
    ax.set_ylabel("Number of Fixed Points")
    ax.scatter(selected_COMs, num_of_XP, s=2, color="#E65100", label="X points")
    ax.scatter(selected_COMs, num_of_OP, s=2, label="O points")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
