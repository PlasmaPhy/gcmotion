from collections import deque
from gcmotion.utils.logger_setup import logger
import matplotlib.pyplot as plt

from gcmotion.utils.bif_values_setup import set_up_bif_plot_values
from gcmotion.configuration.plot_parameters import BifurcationPlotConfig


def _set_up_tpb_base_plot(
    profiles: list,
    COM_plotO: list,
    O_energies_plot: list,
    COM_plotX: list,
    X_energies_plot: list,
    which_COM: str,
    label_energy_units: str,
):
    r"""
    Simple script that sets up the base trapped passing boundary plot
    by checking some parameters.
    """

    if which_COM == "mu":
        x_label_loc = r"$E-{\mu}B_0$" + f"[{label_energy_units}]"
        y_label_loc = r"$\mu$" + f"[{profiles[0].muNU.units}]"
        xO_values = O_energies_plot
        yO_values = COM_plotO
        xX_values = X_energies_plot
        yX_values = COM_plotX
        axis_to_format = "y"

    elif which_COM == "Pzeta":
        x_label_loc = r"$P_{\zeta}$" + f"[{profiles[0].PzetaNU.units}]"
        y_label_loc = f"Energies [{label_energy_units}]"
        xO_values = COM_plotO
        yO_values = O_energies_plot
        xX_values = COM_plotX
        yX_values = X_energies_plot
        axis_to_format = "x"

    else:
        print(
            """\n'which_COM' arguments must be either 'mu' or 'Pzeta'.
            \nABORTING trapped passing boundary plot...\n"""
        )

    return x_label_loc, y_label_loc, xX_values, yX_values, xO_values, yO_values, axis_to_format


def _plot_trapped_passing_boundary(
    profiles: list,
    X_energies: list | deque,
    O_energies: list | deque,
    which_COM: str,
    config: BifurcationPlotConfig,
    ax=None,
):
    r"""Base plotting function. Only draws upon a given axis without showing
    any figures.

    Parameters
    ----------
    profiles : list, deque
        List of profile objects.
    X_energies : deque, list
        The values of the Energies of the X points for each COM value.
    O_energies : deque, list
        The values of the Energies of the O points for each COM value.
    which_COM : str
        String that indicates with respect to which constant of motion (COM) :math:`\mu`
        or :math:`P_{\zeta}` the energies of the fixed points are plotted.
    config : BifurcationPlotConfig
        Dataclass config with several plot parameters. For a full list of all
        available optional parameters, see the dataclass BifurcationPlotConfig at
        gcmotion/configuration/plot_parameters. The defaults values are set there,
        and are overwritten if passed as arguements.
    ax : Axes
        The ax upon which to draw.
    """
    logger.info("\t==> Plotting Base Trapped Passing Boundary...")

    fig_kw = {
        "figsize": config.figsize,
        "dpi": config.dpi,
        "layout": config.layout,
        "facecolor": config.facecolor,
        "sharex": config.sharex,
    }

    fig, ax = plt.subplots(1, 1, **fig_kw)

    # If the selcted COM is mu we need to tilt the energies in the plot by subtracting
    # mu*B0 which is done in set_up_bif_plot_values if tilt_energies = True
    if which_COM == "mu":
        tilted_energies_loc = True
    elif which_COM == "Pzeta":
        tilted_energies_loc = False

    # X O Energies bifurcation plot
    COM_plotX, X_energies_plot = set_up_bif_plot_values(
        profiles=profiles,
        y_values=X_energies,
        which_COM=which_COM,
        tilt_energies=tilted_energies_loc,
        input_energy_units=config.energy_units,
    )
    COM_plotO, O_energies_plot = set_up_bif_plot_values(
        profiles=profiles,
        y_values=O_energies,
        which_COM=which_COM,
        tilt_energies=tilted_energies_loc,
        input_energy_units=config.energy_units,
    )

    x_label_loc, y_label_loc, xX_values, yX_values, xO_values, yO_values, axis_to_format = (
        _set_up_tpb_base_plot(
            profiles=profiles,
            COM_plotO=COM_plotO,
            O_energies_plot=O_energies_plot,
            COM_plotX=COM_plotX,
            X_energies_plot=X_energies_plot,
            which_COM=which_COM,
            label_energy_units=config.energy_units,
        )
    )

    plt.xlabel(x_label_loc)
    ax.set_ylabel(y_label_loc)
    ax.ticklabel_format(style="sci", axis=axis_to_format, scilimits=(0, 0))

    ax.scatter(xX_values, yX_values, s=2, color="#E65100", label="X points")
    ax.scatter(xO_values, yO_values, s=2, label="O points")

    ax.legend()
