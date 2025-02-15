import numpy as np
from collections import deque
from gcmotion.utils.logger_setup import logger

from gcmotion.utils.bif_values_setup import set_up_bif_plot_values


def _thetas_bif_plot(
    profiles: list,
    X_thetas: list | deque,
    O_thetas: list | deque,
    which_COM: str,
    ax=None,
):
    r"""Base plotting function. Only draws upon a given axis without showing
    any figures.

    Parameters
    ----------
    profiles : list, deque
        List of profile objects.
    X_thetas : deque, list
        The values of the :math:`\theta`s of the X points for each COM value.
    O_thetas : deque, list
        The values of the :math:`\theta`s of the O points for each COM value.
    which_COM : str
        String that indicates with respect to which constant of motion (COM) :math:`\mu`
        or :math:`P_{\zeta}` the :math:`\theta`s fixed are plotted.
    ax : Axes
        The ax upon which to draw.
    """
    logger.info("\t==> Plotting Base Thetas Fixed Bifurcation Diagram...")

    # Theta Fixed Bifurcation
    COM_plotX, X_theta_plot = set_up_bif_plot_values(
        profiles=profiles, y_values=X_thetas, which_COM=which_COM
    )
    COM_plotO, O_theta_plot = set_up_bif_plot_values(
        profiles=profiles, y_values=O_thetas, which_COM=which_COM
    )

    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    ax.set_ylabel(r"$\theta_s$ Fixed")
    ax.set_yticks([-np.pi, 0, np.pi])
    ax.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax.set_ylim([-np.pi - 0.5, np.pi + 0.5])
    ax.scatter(COM_plotX, X_theta_plot, s=2, color="#E65100")
    ax.scatter(COM_plotO, O_theta_plot, s=2)
