from scipy.optimize import differential_evolution, fsolve
from collections import namedtuple, deque

from gcmotion.utils.fp_system import system


# Function to locate a single fixed point
def fixed_point(
    method: str,
    initial_condition: list | tuple | deque,
    bounds: list,
    parameters: namedtuple,
    profile: namedtuple,
    known_thetas: bool = False,
    known_theta_value: float = 0,
):

    # System of equations to be solved
    def fixed_point_system(vars):

        if known_thetas:
            psi = vars[0]
            theta = known_theta_value
        else:
            theta, psi = vars

        # We calculate the quantity theta_dot**2+psi_dot**2 or [theta_dot, psi_dot]
        # depending on the method
        system_result = system(
            theta=theta, psi=psi, parameters=parameters, profile=profile, method=method
        )
        return system_result

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
            fixed_point_system, x0=initial_condition, xtol=1e-10, maxfev=5_000, factor=0.1
        )

        psi_solution = result[0] if known_thetas else result[1]
        theta_solution = known_theta_value if known_thetas else result[0]

    else:
        print('method can be either "fsolve" or "differential evolution"')
        return

    return theta_solution, psi_solution
