r""" Function that uses SciPy's :py:func:`differential_evolution` in order to solve 
the algebraic but complicated system of equations :math:`\dot{\theta} = 0 \& \dot{\psi} = 0`

Example
-------

This is how :py:func:`fixed_points` can be called inside the function :py:func:`fixed_points_plot`:

.. code-block:: python

    from gcmotion.scripts.fixed_points import fixed_points as fp

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

    # Calculate fixed points
    _, fixed_points = fp(
        parameters=parameters,
        profile=profile,
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        dist_tol=dist_tol,
        info=info,
    )

    Parameters
    ----------
    parameters : namedtuple
        Namedtuple containing the constants of motion. Currently 
        :math:`\mu, P_{\zeta0}` are used. ATTENTION, magnitude of parameters ONLY must be
        passed, not units as well
    profile : namedtuple
        Dict containing the tokamak configuration objects.
    theta_lim : list, optional
        Provides the limits for the solution search area with regards to the :math:`\theta`
          variable. It will be passed into the "bounds" argument of :py:func:`differential_evolution`. 
    psi_lim : list, optional
        Provides the limits (divided by psi_wall) for the solution search area with regards
        to the :math:`\psi` variable. It will be passed into the "bounds" argument of
        :py:func:`differential_evolution`. 
    dist_tol : float, optional
        Tolerance that determines distinct fixed points. If both :math:`\theta` and
        :math:`\psi` elements of a fixed point are less than :py:data:`dist_tol` apart
        the two fixed points are not considered distinct.
    ic_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\theta` variable 
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution` 
        will be conducted.
    ic_psi_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\psi` variable 
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution` 
        will be conducted.
    info : bool, optional
        Boolean that dictates weather the fixed points and distinct fixed points found 
        will be printed alongside how many where found respectively.

    .. note:: The parameters argument mus contain the parameters in Normalized Units (NU)
    and it must contain their magnitude, NOT the entire Quantity object.

    Returns
    -------

    num_of_dfp, distinct_fixed_points : tuple
        Tuple where the first element is the number of distinct fixed points found and
        the second element is a list containing the distinct points found in the form
        :math:`[\theta_{fixed},\psi_{fixed}]`

"""

import numpy as np
import pint
from itertools import product

from collections import namedtuple
from gcmotion.utils.distinctify import distinctify
from gcmotion.utils.fp_ic_scan import fp_ic_scan as ic_scanner
from gcmotion.utils.single_fixed_point import fixed_point

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def fixed_points(
    parameters: namedtuple,
    profile: namedtuple,
    Q: Quantity,
    method: str = "differential evolution",
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.3],
    dist_tol: float = 1e-3,
    fp_ic_scan_tol: float = 5 * 1e-8,
    ic_theta_grid_density: int = 1000,
    ic_psi_grid_density: int = 1000,
    random_init_cond: bool = False,
    info: bool = False,
    ic_info: bool = False,
    # polish=True,
    # init="sobol",
    # workers=-1,
):

    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

    psi_lim = np.array(psi_lim) * Q("NUpsi_wall")
    psi_lim = psi_lim.to("NUmagnetic_flux").m

    psi_min = psi_lim[0]
    psi_max = psi_lim[1]

    bounds = [(theta_min, theta_max), (0.99 * psi_min, 1.01 * psi_max)]

    if random_init_cond:
        theta_ic = np.linspace(theta_min, theta_max, ic_theta_grid_density)
        psi_ic = np.linspace(psi_min, psi_max, ic_psi_grid_density)

        initial_conditions = [
            [theta_init, psi_init] for theta_init, psi_init in product(theta_ic, psi_ic)
        ]

    else:
        initial_conditions = ic_scanner(
            parameters=parameters,
            profile=profile,
            method=method,
            theta_grid_density=ic_theta_grid_density,
            psi_grid_density=ic_psi_grid_density,
            psi_lim=psi_lim,
            theta_lim=theta_lim,
            tol=fp_ic_scan_tol,
        )

    # Make sure there are initial conditions. If not, use randoom
    if not initial_conditions:
        print("\n\nDID NOT FIND INITIAL CONDITIONS. USING RANDOM\n\n")

        theta_ic = np.linspace(theta_min, theta_max, 3)
        psi_ic = np.linspace(psi_min, psi_max, 70)

        initial_conditions = [
            [theta_init, psi_init] for theta_init, psi_init in product(theta_ic, psi_ic)
        ]

    fixed_points = np.empty((len(initial_conditions), len(initial_conditions[0])))
    fixed_points[:] = np.nan

    idx = 0

    # Run fixed_point() for multiple initial conditions in order to locate
    # multiple fixed points
    for initial_condition in initial_conditions:

        theta_fix, psi_fix = fixed_point(
            method=method,
            initial_condition=initial_condition,
            bounds=bounds,
            parameters=parameters,
            profile=profile,
        )
        fixed_points[idx] = [float(theta_fix), float(psi_fix)]
        idx += 1

    # A lot of the fixed points that were found have identical values-->
    # find out how many distinct fixed points were located
    distinct_fixed_points = distinctify(fixed_points, tol=dist_tol)
    num_of_dfp = distinct_fixed_points.shape[0]

    if ic_info:
        initial_conditions_print = [[float(x), float(y)] for x, y in initial_conditions]
        print(f"\nInitial Conditions:{initial_conditions_print}\n")
        print(f"\nNumber of Initial Conditions: {len(initial_conditions_print)}\n")

    if info and fixed_points.shape[0] <= 30:
        print(f"\nFixed Points: {fixed_points}\n")
        print(f"Number of Fixed Points: {fixed_points.shape[0]}")

    if info:
        print(f"\nDistinct Fixed Points: {distinct_fixed_points}\n")
        print(f"Number of Distinct Fixed Points: {distinct_fixed_points.shape[0]}\n")

    return num_of_dfp, distinct_fixed_points, initial_conditions
