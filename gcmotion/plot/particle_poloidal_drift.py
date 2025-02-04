r"""
Plots the particle's poloidal drift on top of the energy contour plot.
"""

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
    r"""Plots the particle's poloidal drift on top of the energy contour
    plot.

    Parameters
    ----------
    particle: Particle
        The current working particle

    Other Parameters
    ----------------
    thetalim: list, optional
        The :math:`\theta` span. Ignored when plotting in polar projection.
        Defaults to [-π, π].
    psilim: list, optional
        The :math:`\psi` span, relative to :math:`\psi_{wall}`. Defaults to [0,
        1.2]
    projection: {None, "polar"}, optional
        The plot projection. None means cartesian. Defaults to None.
    levels: int, optional
        The energy contour levels. Defaults to 30.
    E_units: str, optional
        The Energy units, as read by pint. Defaults to "keV".
    flux_units: str, optional
        The magnetic flux units, as read by pint. Defaults to "Tesla *
        meter^2".
    grid_density: int, optional
        The contour's grid density. Defaults to 200.
    potential: bool, optional
        Whether or not to add the electric potential in the calculation of
        energy. Defaults to True.
    wall: bool, optional
        Whether or not to shade the area :math:`\psi>\psi_{wall}`. Defaults to
        True.
    cursor: bool, optional
        Whether or not to add cursor data to the plot. This is a mild
        performance hit. Defaults to True.

    """

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
