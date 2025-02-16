import numpy as np
from collections import deque
from gcmotion.utils.logger_setup import logger

from gcmotion.entities.tokamak import Tokamak
from gcmotion.entities.profile import Profile
from gcmotion.utils.bif_values_setup import set_up_bif_plot_values

from gcmotion.configuration.plot_parameters import ParabolasPlotConfig


def _create_profiles_list(profile: Profile, Pzetalim: list | tuple, TPB_density: int):

    tokamak = Tokamak(
        R=profile.R,
        a=profile.a,
        qfactor=profile.qfactor,
        bfield=profile.bfield,
        efield=profile.efield,
    )

    psip_wallNU = profile.psip_wallNU.m

    # Pzetalim is given in PzetaNU/psip_walNU so we need to multiply bu psip_wallNU
    Pzetamin, Pzetamax = Pzetalim  # PzetaNU/psip_wallNU
    Pzetamin *= psip_wallNU  # Pzeta in NUcanmom
    Pzetamax *= psip_wallNU  # Pzeta in NUcanmom

    PzetasNU = np.linspace(Pzetamin, Pzetamax, TPB_density)

    profiles = deque([])

    for PzetaNU in PzetasNU:

        current_profile = Profile(
            tokamak=tokamak,
            species=profile.species,
            mu=profile.muNU,
            Pzeta=profile.Q(PzetaNU, "NUCanonical_momentum"),
        )

        profiles.append(current_profile)

    return list(profiles)


def _plot_parabolas_tpb(
    profile: list,
    X_energies: list | deque,
    O_energies: list | deque,
    x_TPB: np.ndarray,
    ax=None,
    **kwargs,
):
    r"""Base plotting function. Only draws upon a given axis without showing
    any figures.

    Parameters
    ----------
    X_energies : deque, list
        The values of the Energies of the X points for each Pzeta value.
    O_energies : deque, list
        The values of the Energies of the O points for each Pzeta value.
    ax : Axes
        The axes upon which the plot will be realized.

    """
    logger.info("\t==> Plotting Base Parabolas Trapped Passing Boundary...")

    # Unpack parameters
    config = ParabolasPlotConfig()
    for key, value in kwargs.items():
        setattr(config, key, value)

    if not profile.bfield.plain_name == "LAR":

        which_COM = "Pzeta"
        input_energy_units = "NUJoule"
        tilt_energies = False

        profiles = _create_profiles_list(
            profile=profile, Pzetalim=x_TPB, TPB_density=config.TPB_density
        )

        # X O Energies bifurcation plot
        Pzeta_plotX, X_energies_plot = set_up_bif_plot_values(
            profiles=profiles,
            y_values=X_energies,
            which_COM=which_COM,
            tilt_energies=tilt_energies,
            input_energy_units=input_energy_units,
        )
        Pzeta_plotO, O_energies_plot = set_up_bif_plot_values(
            profiles=profiles,
            y_values=O_energies,
            which_COM=which_COM,
            tilt_energies=tilt_energies,
            input_energy_units=input_energy_units,
        )

        # The axes we will be plotting upon have y=E/(μΒ0) and x = Pζ/ψpw
        # but here we have E in NUJoule and Pζ in NUCanmom
        psip_wallNU = profile.psip_wallNU.m
        muNU = profile.muNU.m
        B0NU = profile.bfield.B0.to("NUTesla").m
        muB0NU = muNU * B0NU

        Pzeta_plotO = np.array(Pzeta_plotO) / psip_wallNU  # NUCanmom/ψpw
        Pzeta_plotX = np.array(Pzeta_plotX) / psip_wallNU  # NUCanmom/ψpw

        O_energies_plot = np.array(O_energies_plot) / muB0NU  # NUJoule/(μΒ0)
        X_energies_plot = np.array(X_energies_plot) / muB0NU  # NUJoule/(μΒ0)

        ax.plot(
            Pzeta_plotX,
            X_energies_plot,
            linestyle=config.TPB_X_linestyle,
            color=config.TPB_X_color,
            linewidth=config.linewidth,
            label="X points",
        )
        ax.plot(
            Pzeta_plotO,
            O_energies_plot,
            linestyle=config.TPB_O_linestyle,
            color=config.TPB_O_color,
            linewidth=config.linewidth,
            label="O points",
        )

    else:
        ax.plot(
            x_TPB,
            O_energies,
            linestyle=config.TPB_O_linestyle,
            color=config.TPB_O_color,
            label="TPB_O",
            linewidth=config.linewidth,
        )

        ax.plot(
            x_TPB,
            X_energies,
            linestyle=config.TPB_X_linestyle,
            color=config.TPB_X_color,
            label="TPB_X",
            linewidth=config.linewidth,
        )

    # In case the Pzeta limits are very close to zero
    ax.set_xlim([1.1 * config.Pzetalim[0], 1.1 * abs(config.Pzetalim[1])])
