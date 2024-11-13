r"""
Calculates and plots the fixed points of the system.

Example
-------

This is how :py:func:`fixed_points_plot` can be called inside the function :py:func:`energy_contour`:

.. code-block:: python

    if plot_fixed_points:
        fixed_points_plot(
            cwp,
            theta_lim=[theta_min, theta_max],
            psi_lim=[psi_min, psi_max],
            info=True,
            _internal_call=True,
            canvas=canvas,
            **params,
        )


"""

from gcmotion.scripts.fixed_points import fixed_points as fp
import matplotlib.pyplot as plt
from gcmotion.utils._logger_setup import logger


def fixed_points_plot(
    cwp,
    theta_lim: list,
    psi_lim: list,
    info: bool = False,
    **params,
):
    r"""Draws fixed points plot.

    This method is called internally by ``countour_energy()``
    as well.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    theta_lim : list
        List containing the limits for :math:`\theta` (x- axis limits).
    psi_lim : list
        List containing the limits for :math:`P_{\theta}` (y- axis limits).
    info : bool, optional
        Passed into ``fixed_poits()``. Determines weather to print information
        about the fixed points (number, values). Defaults to ``False``.
    params : dict, optional
        Extra arguements if called for many particles:
    """

    # Unpack params
    _internal_call = params.pop("_internal_call", False)  # POP!
    canvas = params.pop("canvas", None)  # POP!

    if _internal_call:
        logger.disable("gcmotion")
    else:
        logger.enable("gcmotion")

    logger.info(f"Plotting fixed points")

    # Set up fixed points parameters
    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

    P_theta_min = psi_lim[0]
    P_theta_max = psi_lim[1]

    # Set up constants dict
    constants = {"mu": cwp.mu, "mass": cwp.mi, "qi": cwp.qi, "Pzeta0": cwp.Pzeta0}
    # Set up tokamak profile dict
    profile = {
        "qfactor": cwp.qfactor,
        "Bfield": cwp.Bfield,
        "Efield": cwp.Efield,
        "Volts_to_NU": cwp.Volts_to_NU,
    }

    # Calculate fixed points
    _, fixed_points = fp(
        constants,
        profile,
        theta_lim=[theta_min, theta_max],
        psi_lim=[P_theta_min, P_theta_max],
        info=info,
    )

    # Check canvas
    if canvas is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        canvas = (fig, ax)
        logger.debug("\tCreating a new canvas.")
    else:
        fig, ax = canvas
        logger.debug("\tUsing existing canvas.")

    # Extract and plot the [theta, P_theta] fixed points
    thetas_fixed = fixed_points[:, 0]
    P_thetas_fixed = fixed_points[:, 1]

    P_theta_plot = P_thetas_fixed / cwp.psi_wall
    ax.scatter(thetas_fixed, P_theta_plot, marker="x", color="green")
