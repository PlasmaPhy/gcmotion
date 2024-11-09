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
            P_theta_lim=[P_theta_min, P_theta_max],
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
    P_theta_lim : list, optional
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
from scipy.optimize import differential_evolution
from math import sqrt
from gcmotion.utils.distinctify import distinctify


def fixed_points(
    constants: dict,
    profile: dict,
    iterations: int = 50,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    P_theta_lim: list = [0.01, 1.3],
    info: bool = False,
    scaled_P_thetas: bool = False,
    # polish=True,
    # init="sobol",
    # workers=-1,
):

    # Tokamak profile
    qfactor = profile["qfactor"]
    Bfield = profile["Bfield"]
    Efield = profile["Efield"]
    Volts_to_NU = float(profile["Volts_to_NU"])

    # Constants of motion
    # E = float(constants["E"]) # Not used
    mu = float(constants["mu"])
    mi = float(constants["mass"])
    qi = int(constants["qi"])
    Pzeta0 = float(constants["Pzeta0"])  # Not used

    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

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

            # Intermediate values
            phi_der_psip, phi_der_theta = Efield.Phi_der(P_theta)
            phi_der_psip *= Volts_to_NU
            phi_der_theta *= Volts_to_NU
            B_der_psi, B_der_theta = Bfield.B_der(P_theta, theta)
            q = qfactor.q_of_psi(P_theta)
            r = sqrt(2 * P_theta)
            B = Bfield.B(r, theta)
            rho = (Pzeta + qfactor.psip_of_psi(P_theta)) / (Bfield.g)
            par = mu + (qi**2 / mi) * rho**2 * B
            bracket1 = -par * q * B_der_psi + qi * phi_der_psip
            bracket2 = par * B_der_theta + qi * phi_der_theta
            D = Bfield.g * q + Bfield.I

            # Canonical Equations
            theta_dot = (qi / mi) / D * rho * B**2 + Bfield.g / (D * qi) * bracket1
            psi_dot = -Bfield.g * q / (D * qi) * bracket2
            rho_dot = psi_dot / (Bfield.g * q)

            P_theta_dot = psi_dot + rho_dot * Bfield.I

            # -----------------------------------

            # Phi_der_P_theta, _ = Efield.Phi_der(P_theta)
            # Phi_der_P_theta *= Volts_to_NU

            # theta_dot = (
            #     -(1 - sqrt(2 * P_theta) * np.cos(theta))
            #     * (Pzeta + qfactor.psip_of_psi(P_theta)) ** 2
            #     / (sqrt(2 * P_theta))
            #     - mu * np.cos(theta) / (sqrt(2 * P_theta))
            #     + (1 - sqrt(2 * P_theta) * np.cos(theta)) ** 2
            #     * (Pz + qfactor.psip_of_psi(P_theta))
            #     / qfactor.q_of_psi(P_theta)
            #     + Phi_der_P_theta
            # )

            # P_theta_dot = -(1 - sqrt(2 * P_theta) * np.cos(theta)) * (
            #     Pz + qfactor.psip_of_psi(P_theta)
            # ) ** 2 * (sqrt(2 * P_theta) * np.sin(theta)) - mu * (
            #     sqrt(2 * P_theta) * np.sin(theta)
            # )

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
