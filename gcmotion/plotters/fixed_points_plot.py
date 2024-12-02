r"""
Calculates and plots the fixed points of the system.

Example
-------

This is how :py:func:`fixed_points_plot` can be called inside the function :py:func:`energy_contour`:

.. code-block:: python

    if plot_fixed_points:

        params.pop("theta_lim", None)
        params.pop("psi_lim", None)

        fixed_points_plot(
            cwp,
            theta_lim=[-1.01 * np.pi, 1.01 * np.pi],
            psi_lim=psi_lim,
            ic_theta_grid_density=ic_theta_grid_density,
            ic_psi_grid_density=ic_psi_grid_density,
            _internal_call=True,
            dist_tol=dist_tol,
            info=True,
            canvas=canvas,
            **params,
        )


"""

from gcmotion.scripts.fixed_points import fixed_points as fp
from gcmotion.utils.XO_points_classification import XO_points_classification as xoc
import matplotlib.pyplot as plt
import numpy as np

from time import time
from collections import namedtuple
from gcmotion.utils.logger_setup import logger


def fixed_points_plot(
    cwp,
    theta_lim: list,
    psi_lim: list,
    dist_tol: float = 1e-3,
    fp_ic_scan_tol: float = 5 * 1e-8,
    ic_theta_grid_density: int = 1000,
    ic_psi_grid_density: int = 1000,
    random_init_cond: bool = False,
    info: bool = False,
    ic_info: bool = False,
    **params,
):
    r"""Draws fixed points plot.

    This method is called internally by ``energy_contour()``
    as well.

    :meta public:

    Parameters
    ----------
    cwp : :py:class:`~gcmotion.classes.particle.Particle`
        The current working particle.
    theta_lim : list
        List containing the limits for :math:`\theta` (x- axis limits).
    psi_lim : list
        List containing the limits for :math:`\psi` and
        consequently :math:`\psi` (y- axis limits).
    dist_tol : float, optional
        Tolerance that determines distinct fixed points. If both :math:`\theta` and
        :math:`\psi` elements of a fixed point are less than :py:data:`dist_tol` apart
        the two fixed points are not considered distinct.
    ic_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\theta` variable
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution`
        will be conducted. Will be passed to :py:func:`fixed_points`.
    ic_psi_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\psi` variable
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution`
        will be conducted.  Will be passed to :py:func:`fixed_points`.
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

    # Get Tokamak profile
    qfactor = cwp.qfactor
    bfield = cwp.bfield
    efield = cwp.efield

    Profile = namedtuple("Tokamak_Profile", ["qfactor", "bfield", "efield"])
    profile = Profile(
        qfactor=qfactor,
        bfield=bfield,
        efield=efield,
    )

    Parameters = namedtuple("Orbit_Parameters", ["Pzeta0", "mu"])

    # Get Particle Parameters
    parameters = Parameters(
        Pzeta0=cwp.Pzeta0NU.magnitude,
        mu=cwp.muNU.magnitude,
    )

    start = time()
    # Calculate fixed points
    _, fixed_points = fp(
        parameters=parameters,
        profile=profile,
        Q=cwp.Q,
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        dist_tol=dist_tol,
        fp_ic_scan_tol=fp_ic_scan_tol,
        ic_theta_grid_density=ic_theta_grid_density,
        ic_psi_grid_density=ic_psi_grid_density,
        random_init_cond=random_init_cond,
        info=info,
        ic_info=ic_info,
    )
    print(f"\n FIXED POINTS RUN IN {(time() - start):.1f}s\n")

    # CAUTION: The xoc function takes in psis_fixed but returns P_thetas_fixed
    X_points, O_points = xoc(
        unclassified_fixed_points=fixed_points,
        parameters=parameters,
        profile=profile,
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

    # Convert deque to numpy arrays for easy manipulation
    X_thetas, X_P_thetas = zip(*X_points) if X_points else ([], [])
    O_thetas, O_P_thetas = zip(*O_points) if O_points else ([], [])

    X_P_thetas = cwp.Q(X_P_thetas, "NUmagnetic_flux").to("Magnetic_flux")
    O_P_thetas = cwp.Q(O_P_thetas, "NUmagnetic_flux").to("Magnetic_flux")

    ax.set_xticks([-np.pi, 0, np.pi])
    ax.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax.set_xlim([-np.pi, np.pi])
    ax.scatter(X_thetas, X_P_thetas, marker="x", color="#80FF80", s=100)
    ax.scatter(O_thetas, O_P_thetas, marker="o", edgecolor="yellow", facecolors="none", s=100)
