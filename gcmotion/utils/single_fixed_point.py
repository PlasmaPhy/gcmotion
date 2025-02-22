r"""
Script/function that calculates (a single) fixed point of the GC Hamiltonian, based on a
single initial condition.
"""

import numpy as np
from scipy.optimize import differential_evolution, fsolve
from collections import deque
from gcmotion.entities.profile import Profile
from gcmotion.utils.fp_system import system


# Function to locate a single fixed point
def fixed_point(
    method: str,
    initial_condition: list | tuple | deque,
    bounds: list,
    profile: Profile,
    known_thetas: bool = False,
    known_theta_value: float = 0,
) -> tuple[float, float]:
    r"""
    Function that finds a single fixed point of the GC Hamiltonian for a given profile
    (tokamak, particle information). It uses numerical solvers/methods "fsolve" or
    "differential evolution" to locate for which values of the :math:`\theta` and :math:`\psi`
    variables, their time derivatives become zero.


        Parameters
        ----------
        method : str
            Indicates which numerical method (solver) is used. Can be "fsolve" or "differential evolution".
        initial_condition : list | tuple | deque
            Initial condition passed into the numerical solver, in order to begin its search/
            iterative process.
        bounds : list
            Necessary to define a search region for the "differential evolution" numerical
            solver. Not used in "fsolve".
        profile : Profile
            Profile object that contains Tokamak and Particle information.
        known_thetas : bool, optional
            Boolean that indicates weather or not the LAR thetas are used (angles for which we
            know fixed oints occur). Defaults to False
        known_theta_value : float, optional
            Value of the LAR theta used. For this value the numerical solver will search for
            fixed :math:`\psi`. Defaults to 0.

        Returns
        -------
        tuple
            Fixed point (tuple) of the form (:math:`\theta`,:math:`\psi`) that essentially
            represents a solution of the numercial solver. :math:`\psi` is calculated in
            "NUMagnetic_flux" units.

    """

    # System of equations to be solved
    def fixed_point_system(vars):

        if known_thetas:
            psi = vars[0]
            theta = known_theta_value
        else:
            theta, psi = vars

        # We calculate the quantity theta_dot**2+psi_dot**2 or [theta_dot, psi_dot]
        # depending on the method
        system_result = system(theta=theta, psi=psi, profile=profile, method=method)

        if known_thetas and method == "fsolve":
            return system_result[0] ** 2 + system_result[1] ** 2

        return system_result

    # Jacobian of equations
    def fixed_point_system_jacobian(vars):

        if known_thetas:
            psi = vars[0]
            theta = known_theta_value
        else:
            theta, psi = vars

        # We calculate the Jacobian
        jacobian = _jacobian(theta=theta, psi=psi, profile=profile)

        return jacobian

    if method == "differential evolution":

        result = differential_evolution(
            fixed_point_system,
            bounds,
            x0=initial_condition,
            tol=1e-9,
            atol=1e-15,
            maxiter=15_000,
            popsize=20,  # 15,
            mutation=(0.5, 1),  # (0.7, 1.5),
            recombination=0.7,  # 0.8,
            strategy="best1bin",  # "best2bin",
        )

        psi_solution = result.x[0] if known_thetas else result.x[1]
        theta_solution = known_theta_value if known_thetas else result.x[0]

    elif method == "fsolve":

        result = fsolve(
            fixed_point_system,
            x0=initial_condition,
            xtol=1e-10,
            fprime=fixed_point_system_jacobian,
            maxfev=2_00,
            factor=0.1,
        )

        psi_solution = result[0] if known_thetas else result[1]
        theta_solution = known_theta_value if known_thetas else result[0]

    else:
        print('method can be either "fsolve" or "differential evolution"')
        return

    return theta_solution, psi_solution


def _jacobian(theta: float, psi: float, profile: Profile) -> np.ndarray:

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
    mu = profile.muNU.m
    Pzeta0 = profile.PzetaNU.m
    Pzeta = Pzeta0

    psi = max(
        psi, profile.psi_wallNU.m / 100
    )  # Should become small for Pzetas close to 0 because psi-->0

    # Object methods calls
    q, q_der_psi = solverqNU(psi)
    psi_p = psipNU(psi)
    b, b_der, currents, currents_der = solverbNU(psi, theta)
    phi_der_psi, phi_der_theta, d2Phi_dpsi2 = solverPhiderNU(psi, theta)

    # Unpack
    b_der_psi, b_der_theta, d2b_dpsi2, d2b_dtheta2, d2b_dpsidtheta = b_der
    i, g = currents
    i_der, g_der, d2i_dpsi2, d2g_dpsi2 = currents_der  # CAUTION: Derivatives with respect to psi

    # Intermediate values
    rho = (Pzeta + psi_p) / g
    par = mu + rho**2 * b
    bracket1 = par * b_der_psi + phi_der_psi
    bracket2 = par * b_der_theta + phi_der_theta
    D = g * q + i + rho * q * (g * i_der - i * g_der)

    gq_D = g * q / D

    D_der_psi = (
        (g_der * q + g * q_der_psi)
        + i_der
        + rho * ((q_der_psi * (g * i_der - i * g_der)) + (q * (g * d2i_dpsi2 - i * d2g_dpsi2)))
    )

    gq_D_der_psi = (1 / D) * (g_der * q + g * q_der_psi - g * q * D_der_psi / D)

    # Jacobian Element Terms

    J11a = 2 * rho * b * b_der_theta / D * (1 - rho * q * g_der)
    J11b = gq_D * (rho**2 * b_der_theta * b_der_psi + par * d2b_dpsidtheta)

    J12a = (rho / D) * (2 * b * b_der_psi - b**2 * D_der_psi / D) * (1 - rho * q * g_der)
    J12b = -(rho**2 * b**2 / D) * (q_der_psi * g_der + q * d2g_dpsi2)
    J12c = gq_D_der_psi * bracket1
    J12d = gq_D * (rho**2 * b_der_psi**2 + par * d2b_dpsi2 + d2Phi_dpsi2)

    J21a = -gq_D * (rho**2 * b_der_theta**2 + par * d2b_dtheta2)

    J22a = -gq_D_der_psi * bracket2
    J22b = -gq_D * (rho**2 * b_der_psi * b_der_theta + par * d2b_dpsidtheta)

    # Jacobian Element

    J11 = J11a + J11b
    J12 = J12a + J12b + J12c + J12d
    J21 = J21a
    J22 = J22a + J22b

    return np.array([[J11, J12], [J21, J22]])
