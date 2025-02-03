import matplotlib.pyplot as plt

from gcmotion.entities.particle import Particle
from gcmotion.configuration.plot_parameters import ParticlePoloidalDrift
from gcmotion.plot._base._base_profile_Energy_contour import (
    _base_profile_Energy_contour,
)
from gcmotion.plot._base._base_contour_colorbar import _base_contour_colorbar
from gcmotion.plot._base._base_particle_poloidal_drift import (
    _base_particle_poloidal_drift,
)


def particle_poloidal_drift(particle: Particle, **kwargs):

    # Unpack parameters
    config = ParticlePoloidalDrift()
    for key, value in kwargs.items():
        setattr(config, key, value)

    # Create figure
    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
        "facecolor": config.facecolor,
    }
    fig = plt.figure(**fig_kw)
    driftax = fig.add_subplot(projection=config.projection)

    # ==============
    # Energy Contour
    # ==============
    # Draw the contour and get the contour object
    Contour = _base_profile_Energy_contour(
        profile=particle.profile, ax=driftax, **kwargs
    )

    # ==============
    # Poloidal Drift
    # ==============
    _base_particle_poloidal_drift(particle=particle, ax=driftax, **kwargs)

    # ========
    # Colorbar
    # ========
    # 'cax=None' creates a new axes for the colorbar, by stealing space from
    # the `ax=contour`
    cbar = fig.colorbar(Contour, cax=None, ax=driftax)

    # Now that the colorbar is created, pass its "ax" to be customized
    _base_contour_colorbar(ax=cbar.ax, contour=Contour, numticks=10)

    cbar.ax.set_title(
        label=f"Energy [{config.E_units}]", size=config.cbarlabelsize
    )

    driftax.set_ylim(
        particle.Q(config.psilim, "NUpsi_wall").to(config.flux_units)
    )

    show = kwargs.get("show", True)
    if show:
        plt.show()
