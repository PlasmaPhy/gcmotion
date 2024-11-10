r"""
Plots the q-factor, magnetic and electric profile of the tokamak.

The zoom option zooms every plot in the :math:`\psi`-axis, 
relative to :math:`\psi_{wall}`. 

Example
-------

.. code-block:: python

    gcm.tokamak_profile(cwp, zoom=[0.7, 1.1])

.. rubric:: Function:
    :heading-level: 4

"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.logger_setup import logger

from gcmotion.configuration.plot_parameters import tokamak_profile as config


def tokamak_profile(cwp, zoom: list = [0, 1.1]):
    r"""Plots the electric field, potential, and q factor,
    with respect to :math:`\psi/\psi_{wall}`.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    zoom : list, optional
        Zoom to specific area in the x-axis of the electric field
        and potential plots. Defaults to [0, 1.1].
    """
    logger.info("Plotting tokamak profile...")

    # Get all needed attributes first
    psi_wall = cwp.psi_wall
    qfactor = cwp.qfactor
    bfield = cwp.bfield
    efield = cwp.efield

    fig = plt.figure(figsize=(12, 12))
    fig.subplots_adjust(hspace=0.5)
    ax_phi = fig.add_subplot(321)
    ax_E = fig.add_subplot(322)
    ax_q1 = fig.add_subplot(323)
    ax_q2 = fig.add_subplot(324)

    zoom = np.array(zoom) * psi_wall
    psis = np.linspace(zoom[0], zoom[1], 1000)

    def plot_electric():
        """Plots the electric field profile in subplots 321 and 322."""
        logger.debug("\tPlotting electric field profile...")
        nonlocal psis
        Er = efield.Er(psis.magnitude)
        Phi = efield.PhiNU(psis.magnitude)
        if np.abs(Er).max() > 1000:  # If in kV
            Er /= 1000
            Phi /= 1000
            E_ylabel = "$E_r$ [kV/m]"
            Phi_ylabel = "$Φ_r$ [kV]"
        else:  # If in V
            E_ylabel = "$E_r$ [V/m]"
            Phi_ylabel = "$Φ_r$ [V]"

        # Radial E field
        ax_phi.plot(psis / psi_wall, Er, color="b", linewidth=1.5)
        ax_phi.plot([1, 1], [Er.min(), Er.max()], color="r", linewidth=1.5)
        ax_phi.set_xlabel(r"$\psi/\psi_{wall}$")
        ax_phi.set_ylabel(E_ylabel)
        ax_phi.set_title("Radial electric field [kV/m]", c="b")

        # Electric Potential
        ax_E.plot(psis / psi_wall, Phi, color="b", linewidth=1.5)
        ax_E.plot([1, 1], [Phi.min(), Phi.max()], color="r", linewidth=1.5)
        ax_E.set_xlabel(r"$\psi/\psi_{wall}$")
        ax_E.set_ylabel(Phi_ylabel)
        ax_E.set_title("Electric Potential [kV]", c="b")

        # ax_phi.set_xlim(zoom)
        # ax_E.set_xlim(zoom)

        logger.debug("\t-> Electric field profile successfully plotted.")

    def plot_q():
        """Plots the q factor profile in subplots 323 and 324."""
        logger.debug("\tPlotting q factor profile...")
        nonlocal psis
        # q(ψ)
        y1 = qfactor.q_of_psi(psis.magnitude)
        if type(y1) is int:  # if q = Unity
            y1 *= np.ones(psis.shape)
        ax_q1.plot(psis / psi_wall, y1, color="b", linewidth=1.5)
        ax_q1.plot([1, 1], [y1.min(), y1.max()], color="r", linewidth=1.5)
        ax_q1.set_xlabel(r"$\psi/\psi_{wall}$")
        ax_q1.set_ylabel(r"$q(\psi)$")
        ax_q1.set_title(r"$\text{q factor }q(\psi)$", c="b")

        # ψ_π(ψ)
        y2 = qfactor.psipNU(psis.magnitude)
        ax_q2.plot(psis / psi_wall, y2, color="b", linewidth=1.5)
        ax_q2.plot([1, 1], [y2.min(), y2.max()], color="r", linewidth=1.5)
        ax_q2.set_xlabel(r"$\psi/\psi_{wall}$")
        ax_q2.set_ylabel(r"$\psi_p(\psi)$")
        ax_q2.set_title(r"$\psi_p(\psi)$", c="b")

        logger.debug("\t-> Q factor profile successfully plotted.")

    def plot_magnetic():
        """Plots the magnetic field profile in a single bottom subplot."""
        logger.debug("\tPlotting electric field profile...")
        nonlocal psis
        ax_B = plt.subplot(212, projection="polar")
        box = ax_B.get_position()
        box.y0 = box.y0 - 0.1
        box.y1 = box.y1 - 0.1
        ax_B.set_position(box)
        ax_B.set_title(
            "LAR Magnetic Field Profile", c="b", loc="center", pad=5
        )
        ax_B.set_rlabel_position(30)
        low_ylim = max(0, zoom[0].magnitude)
        yticks = (
            np.array([low_ylim, zoom[1].magnitude, psi_wall.magnitude])
            * zoom.units
            / psi_wall
        )
        ax_B.yaxis.set_units(yticks.units)
        ax_B.set_yticks(yticks)

        thetas = np.linspace(0, 2 * np.pi, len(psis))
        psi, theta = np.meshgrid(psis.magnitude, thetas)
        B = bfield.b(psi, theta)
        levels = config["contour_params"]["levels"]
        cmap = config["contour_params"]["cmap"]
        ax_B.contourf(
            theta, psi / psi_wall.magnitude, B, levels=levels, cmap=cmap
        )

        logger.debug("\t-> Magnetic field profile successfully plotted.")

    plot_electric()
    plot_q()
    plot_magnetic()

    # fig.set_tight_layout(True)
    # plt.tight_layout()
    plt.ion()
    plt.show(block=True)

    logger.info("--> Tokamak profile successfully plotted.")
