import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import LogLocator
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from gcmotion.entities.profile import Profile
from gcmotion.configuration.scripts_configuration import ContourFreqConfig

tau = 2 * np.pi


def _ptheta_energy_contour(
    profile: Profile, psilim: tuple, config: ContourFreqConfig
):
    r"""Creates the needed θ-Pθ contour to calculate ωθ frequencies.

    Parameters
    ----------
    profile: Profile
        The Profile object upon which the analysis is performed
    psilim: tuple
        The grid ψ limit to calculate the energy on. Note that ψ is only used
        to calculate the energy; Afterwards, the grid contains the
        corresponding Pθ values.
    config: ContourFreqConfig
        Dataclass containing configuration options.

    """

    thetalim = profile.Q((-tau, tau), "radians")
    psilim = profile.Q(psilim, "psi_wall").to("NUMagnetic_flux")
    psi_grid, theta_grid = np.meshgrid(
        np.linspace(psilim[0], psilim[1], config.grid_density),
        np.linspace(thetalim[0], thetalim[1], config.grid_density),
    )

    energy_grid = profile.findEnergy(
        psi=psi_grid,
        theta=theta_grid.m,
        units="NUJoule",
        potential=config.potential,
    )
    ptheta_grid = profile.findPtheta(
        psi=psi_grid,
        units="NUCanonical_momentum",
    )

    data = {
        "theta": theta_grid.m,
        "Ptheta": ptheta_grid.m,
        "Energy": energy_grid.m,
    }

    locator = LogLocator(base=config.log_base, numticks=config.levels)
    locator.MAXTICKS = 10000

    kw = {
        "locator": locator,
    }

    fig, ax = plt.subplots()
    CPtheta = ax.contour("theta", "Ptheta", "Energy", data=data, **kw)

    # Store the data to access them later
    CPtheta.data = data
    CPtheta.ylim = (ptheta_grid.m[0][0], ptheta_grid.m[0][-1])

    if config.show_master_contour:
        _plot_master_contour(ax)
    else:
        plt.close()

    return CPtheta


def _plot_master_contour(fig: Figure, ax: Axes):
    r"""Plot the master contour."""

    fig.set_constrained_layout(True)
    ax.set_title(r"$Master \theta-P_\theta Contour$")
    ax.set_xlabel(r"$\theta [radians]$")
    ax.set_ylabel(r"$P_\theta [NU]$")

    plt.show()
