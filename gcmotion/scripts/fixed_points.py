r""" Function that uses SciPy's :py:func:`differential_evolution` in order to solve 
the algebraic but complicated system of equations :math:`\dot{\theta} = 0 \& \dot{P_{\theta}} = 0`

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
        to the :math:`P_{\theta}` variable. It will be passed into the "bounds" argument of
        :py:func:`differential_evolution`. 
    dist_tol : float, optional
        Tolerance that determines distinct fixed points. If both :math:`P_{\theta}` and
        :math:`P_{\theta}` elements of a fixed point are less than :py:data:`dist_tol` apart
        the two fixed points are not considered distinct.
    ic_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\theta` variable 
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution` 
        will be conducted.
    ic_P_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`P_{\theta}` variable 
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
        :math:`[\theta_{fixed},P_{\theta_{fixed}}]`

"""

import numpy as np
import pint
from scipy.optimize import differential_evolution, fsolve

from collections import namedtuple
from gcmotion.utils.distinctify import distinctify
from gcmotion.utils.energy_Ptheta import energy_Ptheta
from gcmotion.scripts.fp_ic_scan import fp_ic_scan as ic_scanner

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def fixed_points(
    parameters: namedtuple,
    profile: namedtuple,
    Q: Quantity,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.3],
    dist_tol: float = 1e-3,
    ic_theta_grid_density: int = 800,
    ic_P_theta_grid_density: int = 800,
    info: bool = False,
    # polish=True,
    # init="sobol",
    # workers=-1,
):

    # Tokamak profile
    qfactor = profile.qfactor
    bfield = profile.bfield
    efield = profile.efield

    # Define quantites for the solver for clarity
    solverqNU = qfactor.solverqNU
    psipNU = qfactor.psipNU
    solverbNU = bfield.solverbNU
    solverPhiderNU = efield.solverPhiderNU

    # Parameters ()
    mu = parameters.mu
    Pzeta0 = parameters.Pzeta0

    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

    psi_lim = np.array(psi_lim) * Q("NUpsi_wall")
    psi_lim = psi_lim.to("NUmagnetic_flux").m

    _, P_theta_lim = energy_Ptheta(
        psi=psi_lim, theta=0, mu=mu, Pzeta=Pzeta0, profile=profile, contour_Phi=True
    )

    P_theta_min = P_theta_lim[0]
    P_theta_max = P_theta_lim[1]

    bounds = [(theta_min, theta_max), (0.99 * P_theta_min, 1.01 * P_theta_max)]

    # Function to locate a single fixed point
    def fixed_point(initial_condition=None):

        Pzeta = Pzeta0

        # System of equations to be solved
        def system(vars):

            theta, P_theta = vars
            P_theta = max(P_theta, P_theta_min)

            # Object methods calls
            q = solverqNU(P_theta)
            psi_p = psipNU(P_theta)
            b, b_der, currents, currents_der = solverbNU(P_theta, theta)
            phi_der_psi, phi_der_theta = solverPhiderNU(P_theta, theta)

            # Unpack
            b_der_psi, b_der_theta = b_der
            i, g = currents
            i_der, g_der = currents_der
            # Multiply current derivatives by q to get their derivatives with
            # respect to psi instead of psip
            i_der, g_der = q * i_der, q * g_der

            # Intermediate values
            rho = (Pzeta + psi_p) / g
            par = mu + rho**2 * b
            bracket1 = -par * b_der_psi + phi_der_psi
            bracket2 = par * b_der_theta + phi_der_theta
            D = g * q + i + rho * (g * i_der - i * g_der)

            # Canonical Equations
            theta_dot = (1 - rho * g_der) / D * rho * b**2 + q * g / D * bracket1
            psi_dot = -q * g / D * bracket2
            rho_dot = -(1 - rho * g_der) / D * bracket2

            # Canonical Momentum
            P_theta_dot = psi_dot + rho_dot * i

            return theta_dot**2 + P_theta_dot**2

        # Use differential evolution to find solutions
        result = differential_evolution(
            system,
            bounds,
            x0=initial_condition,
            tol=1e-7,
            atol=1e-15,
            maxiter=20_000,
            popsize=15,
            mutation=(0.5, 1),
            recombination=0.7,
            strategy="best1bin",
        )
        theta_solution, P_theta_solution = result.x

        return theta_solution, P_theta_solution

    initial_conditions = ic_scanner(
        parameters=parameters,
        profile=profile,
        theta_grid_density=ic_theta_grid_density,
        P_theta_grid_density=ic_P_theta_grid_density,
        psi_lim=psi_lim,
        theta_lim=theta_lim,
        tol=1e-6,
        info=info,
    )

    fixed_points = np.empty((len(initial_conditions), len(initial_conditions[0])))
    fixed_points[:] = np.nan

    idx = 0

    # Run fixed_point() for multiple initial conditions in order to locate
    # multiple fixed points
    for initial_condition in initial_conditions:

        theta_fix, P_theta_fix = fixed_point(initial_condition=initial_condition)
        fixed_points[idx] = [float(theta_fix), float(P_theta_fix)]
        idx += 1

    # A lot of the fixed points that were found have identical values-->
    # find out how many distinct fixed points were located
    distinct_fixed_points = distinctify(fixed_points, tol=dist_tol)
    num_of_dfp = distinct_fixed_points.shape[0]

    if info:

        print(f"\nFixed Points: {fixed_points}\n")
        print(f"Number of Fixed Points: {fixed_points.shape[0]}")

        print(f"\nDistinct Fixed Points: {distinct_fixed_points}\n")
        print(f"Number of Distinct Fixed Points: {distinct_fixed_points.shape[0]}\n")

    return num_of_dfp, distinct_fixed_points
