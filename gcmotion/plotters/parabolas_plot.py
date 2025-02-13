import numpy as np
from gcmotion.utils.logger_setup import logger
import matplotlib.pyplot as plt

from gcmotion.configuration.parabolas_parameters import ParabolasPlotConfig
from gcmotion.entities.profile import Profile
from gcmotion.scripts.parabolas import calc_parabolas


def parabolas_diagram(profile: Profile, ax=None, **kwargs):

    # Unpack parameters
    config = ParabolasPlotConfig()
    for key, value in kwargs.items():
        setattr(config, key, value)

    # Create figure if ax is not given
    if ax is None:
        fig_kw = {
            "figsize": config.figsize,
            "dpi": config.dpi,
            "layout": config.layout,
            "facecolor": config.facecolor,
        }

        fig, par_ax = plt.subplots(1, 1, **fig_kw)
    else:
        par_ax = ax

    # Calculate parabolas values
    x, y_R, y_L, y_MA, LAR_tpb_O, LAR_tpb_X = calc_parabolas(
        Pzetalim=config.Pzetalim, profile=profile, Pzeta_density=config.Pzeta_density
    )

    # Plot right left boundar, magnetic axis
    par_ax.plot(x, y_R, linestyle="dashed", color="orange", label="RW", linewidth=config.linewidth)
    par_ax.plot(x, y_L, linestyle="dashdot", color="orange", label="LW", linewidth=config.linewidth)
    par_ax.plot(x, y_MA, linestyle="dotted", color="orange", label="MA", linewidth=config.linewidth)

    # If asked and if LAR plot LAR trapped-passing boundary
    if config.LAR_TPB and profile.bfield.plain_name == "LAR":
        par_ax.plot(x, LAR_tpb_O, linestyle="solid", label="TPB_LAR_O", linewidth=config.linewidth)
        par_ax.plot(
            x,
            LAR_tpb_X,
            linestyle="solid",
            color="#E65100",
            label="TPB_LAR_X",
            linewidth=config.linewidth,
        )

    par_ax.set_xlabel(r"$\dfrac{P_\zeta}{\psi_{p_w}}$", fontsize=config.xlabel_fontsize)
    par_ax.set_ylabel(
        r"$\dfrac{E}{\mu B_0}$", rotation=config.ylabel_rotation, fontsize=config.ylabel_fontsize
    )
    par_ax.set_title(config.title)
    par_ax.set_ylim(config.enlim)
    if config.legend:
        par_ax.legend()

    plt.show()
    plt.ion()
