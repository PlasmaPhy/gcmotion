import numpy as np
from collections import deque
from gcmotion.utils.logger_setup import logger

from gcmotion.utils.bif_values_setup import set_up_bif_plot_values


def _P_thetas_bif_plot(
    profiles: list,
    X_P_thetas: list | deque,
    O_P_thetas: list | deque,
    which_COM: str,
    ax=None,
):
    r"""Base plotting function. Only draws upon a given axis without showing
    any figures.

    Parameters
    ----------
    profiles : list, deque
        List of profile objects.
    X_P_thetas : deque, list
        The values of the :math:`P_{\theta}`s of the X points for each COM value.
    O_P_thetas : deque, list
        The values of the :math:`P_{\theta}`s of the O points for each COM value.
    which_COM : str
        String that indicates with respect to which constant of motion (COM) :math:`\mu`
        or :math:`P_{\zeta}` the :math:`P_{\theta}`s fixed are plotted.
    ax : Axes
        The ax upon which to draw.
    """
    logger.info("\t==> Plotting Base P_thetas Fixed Bifurcation Diagram...")

    # Theta Fixed Bifurcation
    COM_plotX, X_P_theta_plot = set_up_bif_plot_values(
        profiles=profiles, y_values=X_P_thetas, which_COM=which_COM
    )
    COM_plotO, O_P_theta_plot = set_up_bif_plot_values(
        profiles=profiles, y_values=O_P_thetas, which_COM=which_COM
    )

    # Set the upper limit of the y axis properly
    # Combine the two lists for comparison
    psi_wallNU = profiles[0].psi_wall.to("NUMagnetic_flux")
    P_theta_wallNU = profiles[0].findPtheta(psi=psi_wallNU, units="NUCanonical_momentum").m

    combined_P_theta_plot = X_P_theta_plot + O_P_theta_plot
    ul = max(combined_P_theta_plot) * 1.05

    if np.isclose(max(combined_P_theta_plot), P_theta_wallNU):
        ul = P_theta_wallNU * 1.05

    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    ax.set_ylabel(r"$P_{\theta_s}$ Fixed")
    ax.set_ylim(0, ul)
    ax.scatter(COM_plotX, X_P_theta_plot, s=2, color="#E65100", label="X points")
    ax.scatter(COM_plotO, O_P_theta_plot, s=2, label="O points")
    ax.axhline(y=P_theta_wallNU, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.legend(loc="lower left")
