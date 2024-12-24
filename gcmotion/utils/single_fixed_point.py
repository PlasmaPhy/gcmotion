from scipy.optimize import differential_evolution, fsolve
from collections import namedtuple

from gcmotion.utils.fp_system import system


# Function to locate a single fixed point
def fixed_point(initial_condition, bounds: list, parameters: namedtuple, profile: namedtuple):

    # System of equations to be solved
    def fixed_point_system(vars):

        theta, psi = vars

        # We calculate the quantity theta_dot**2+psi_dot**2
        theta_dot_sq_psi_dot_sq_sum = system(
            theta=theta, psi=psi, parameters=parameters, profile=profile
        )

        return theta_dot_sq_psi_dot_sq_sum

    # Use differential evolution to find solutions
    # result = differential_evolution(
    #     fixed_point_system,
    #     bounds,
    #     x0=initial_condition,
    #     tol=1e-9,
    #     atol=1e-15,
    #     maxiter=15_000,
    #     popsize=20,  # 15,
    #     mutation=(0.5, 1),  # (0.7, 1.5),
    #     recombination=0.7,  # 0.8,
    #     strategy="best1bin",  # "best2bin",
    # )

    result = fsolve(fixed_point_system, x0=initial_condition, xtol=1e-10, maxfev=5_000, factor=0.1)

    theta_solution, psi_solution = result

    return theta_solution, psi_solution
