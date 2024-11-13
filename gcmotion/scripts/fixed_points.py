r""" Function that uses SciPy's :py:func:`differential_evolution` in order to solve 
the algebraic but complicated system of equations :math:`\dot{\theta} = 0 \& \dot{P_{\theta}} = 0`

Example
-------

This is how :py:func:`fixed_points` can be called inside the function :py:func:`fixed_points_plot`:

.. code-block:: python

    from gcmotion.scripts.fixed_points import fixed_points as fp

        constants = {"mu": cwp.mu, "mass": cwp.mi, "qi": cwp.qi, "Pzeta0": cwp.Pzeta0}

        profile = {
            "qfactor": cwp.qfactor,
            "Bfield": cwp.Bfield,
            "Efield": cwp.Efield,
            "Volts_to_NU": cwp.Volts_to_NU,
        }

        _, fixed_points = fp(
            constants,
            profile,
            theta_lim=[theta_min, theta_max],
            psi_lim=[P_theta_min, P_theta_max],
            info=info,
        )




    Parameters
    ----------
    constants : dict
        Dict containing the constants of motion. Currently :math:`\mu, m_i, q_i, P_{\zeta0}`
        are used.
    profile : dict
        Dict containing the tokamak configuration objects.
    iterations : int, optional
        Integer that essentially dictates the number of initial conditions that will be fed into
        differential_evolution, thus dictating the number of iterations the
        afformentioned function will be executed.
    theta_lim : list, optional
        Provides the limits for the solution search area with regards to the :math:`\theta`
          variable. It will be passed into the "bounds" argument of :py:func:`differential_evolution`. 
    psi_lim : list, optional
        Provides the limits (divided by psi_wall) for the solution search area with regards
        to the :math:`P_{\theta}` variable. It will be passed into the "bounds" argument of
        :py:func:`differential_evolution`. 
    info : bool, optional
        Boolean that dictates weather the fixed points and distinct fixed points found 
        will be printed alongside how many where found respectively.
    scaled_P_thetas : bool, optional
        Boolean that dictates weather the initial conditions regarding the variable
        :math:`P_{\theta}` (initial conditions for fixed points search) will be scaled, so
        as to increase the density of initial :math:`P_{\theta}`'s near the lower boundary
        of the search area.

    Returns
    -------

    num_of_dfp, distinct_fixed_points : tuple
        Tuple where the first element is the number of distinct fixed points found and
        the second element is a list containing the distinct points found in the form
        :math:`[\theta_{fixed},P_{\theta_{fixed}}]`

"""

import numpy as np
import pint
from scipy.optimize import differential_evolution

from collections import namedtuple
from gcmotion.utils.distinctify import distinctify
from gcmotion.utils.energy_Ptheta import energy_Ptheta

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def fixed_points(
    parameters: namedtuple,
    profile: namedtuple,
    Q: Quantity,
    iterations: int = 50,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.3],
    info: bool = False,
    scaled_P_thetas: bool = False,
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

    P_theta_lim = energy_Ptheta(
        psi=psi_lim, theta=0, mu=mu, Pzeta=Pzeta0, profile=profile, contour_Phi=True
    )

    P_theta_min = P_theta_lim[0]
    P_theta_max = P_theta_lim[1]

    bounds = [(theta_min, theta_max), (P_theta_min, P_theta_max)]

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
            maxiter=10000,
            popsize=15,
            mutation=(0.5, 1),
            recombination=0.7,
            strategy="best1bin",
        )
        theta_solution, P_theta_solution = result.x

        return theta_solution, P_theta_solution

    fixed_points = np.empty((3 * iterations, 2))
    fixed_points[:] = np.nan

    # Define the limits of the initial conditions to be certainly inside the
    # limits of the "bounds" argument you have defined previously
    if theta_min <= 0:
        theta_min_init = theta_min * 0.999
    else:
        theta_min_init = theta_min * 1.001

    theta_max_init = np.pi * 0.999

    P_theta_min_init = P_theta_min * 1.001
    P_theta_max_init = P_theta_max * 0.999

    thetas_init = np.linspace(theta_min_init, theta_max_init, 3)
    P_thetas_init = np.linspace(P_theta_min_init, P_theta_max_init, iterations)

    if scaled_P_thetas:
        P_thetas_init_mod = np.sin(np.pi * P_thetas_init) ** 4  # (1 - np.cos(P_thetas_init)) / 2
        P_thetas_init = P_theta_min_init + (P_theta_max_init - P_theta_min_init) * P_thetas_init_mod

    idx = 0

    # Run fixed_point() for multiple initial conditions in order to locate
    # multiple fixed points
    for theta_init in thetas_init:
        for P_theta_init in P_thetas_init:

            initial_condition = [theta_init, P_theta_init]

            theta_fix, P_theta_fix = fixed_point(initial_condition=initial_condition)
            fixed_points[idx] = [float(theta_fix), float(P_theta_fix)]
            idx += 1

    # A lot of the fixed points that were found have identical values-->
    # find out how many distinct fixed points were located
    distinct_fixed_points = distinctify(fixed_points)
    num_of_dfp = distinct_fixed_points.shape[0]

    if info:

        print(f"\nFixed Points: {fixed_points}\n")
        print(f"Number of Fixed Points: {fixed_points.shape[0]}")

        print(f"\nDistinct Fixed Points: {distinct_fixed_points}\n")
        print(f"Number of Distinct Fixed Points: {distinct_fixed_points.shape[0]}\n")

    return num_of_dfp, distinct_fixed_points
