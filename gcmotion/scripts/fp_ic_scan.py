import numpy as np
import pint

from gcmotion.utils.energy_Ptheta import energy_Ptheta
from collections import namedtuple, deque

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def fp_ic_scan(
    parameters: namedtuple,
    profile: namedtuple,
    theta_grid_density: int = 800,
    P_theta_grid_density: int = 800,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.8],
    tol: float = 1e-6,
    info: bool = False,
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

    # WE SKIP THIS STEP BECAUSE IT HAPPENS IN FIXED POINTS ALREADY SO psi_lim PASSED HERE
    # WILL ALREADY BE PASSED THROUGH THESE 2 LINES OF CODE
    # psi_lim = np.array(psi_lim) * Q("NUpsi_wall")
    # psi_lim = psi_lim.to("NUmagnetic_flux").m

    _, P_theta_lim = energy_Ptheta(
        psi=psi_lim, theta=0, mu=mu, Pzeta=Pzeta0, profile=profile, contour_Phi=True
    )

    P_theta_min = 0.99 * P_theta_lim[0]
    P_theta_max = 1.01 * P_theta_lim[1]

    # Create the grid
    theta_grid, P_thetaNU_grid = np.meshgrid(
        np.linspace(theta_min, theta_max, theta_grid_density),
        np.linspace(P_theta_min, P_theta_max, P_theta_grid_density),
    )

    grid_shape = theta_grid.shape

    def system(theta, P_theta):

        if P_theta < 0:
            print("WARNING: You have given values resulting in P_theta < 0")
            return

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
        rho = (Pzeta0 + psi_p) / g
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

    # List to store the coordinates of points where the quantity becomes zero
    minima = deque([])

    for i in range(1, grid_shape[0] - 1):  # Exclude first and last row
        for j in range(1, grid_shape[1] - 1):  # Exclude first and last column
            # Evaluate the system at the current grid point

            center = system(theta_grid[i, j], P_thetaNU_grid[i, j])

            # If the value of theta_dot**2+P_theta_dot**2 is not sufficiently
            # small to be considered 0 we do not have a fixed point candidate
            if center >= tol:
                continue

            # Evaluate the 8 neighbors
            neighbors = [
                system(theta_grid[i - 1, j - 1], P_thetaNU_grid[i - 1, j - 1]),
                system(theta_grid[i - 1, j], P_thetaNU_grid[i - 1, j]),
                system(theta_grid[i - 1, j + 1], P_thetaNU_grid[i - 1, j + 1]),
                system(theta_grid[i, j - 1], P_thetaNU_grid[i, j - 1]),
                system(theta_grid[i, j + 1], P_thetaNU_grid[i, j + 1]),
                system(theta_grid[i + 1, j - 1], P_thetaNU_grid[i + 1, j - 1]),
                system(theta_grid[i + 1, j], P_thetaNU_grid[i + 1, j]),
                system(theta_grid[i + 1, j + 1], P_thetaNU_grid[i + 1, j + 1]),
            ]

            # Check if the center value is smaller than all neighbors
            if all(center < n for n in neighbors):
                minima.append((theta_grid[i, j], P_thetaNU_grid[i, j]))

    initial_conditions = minima.copy()

    if info:
        initial_conditions_pr = [[float(x), float(y)] for x, y in initial_conditions]
        print(f"\nInitial Conditions:{initial_conditions_pr}\n")
        print(f"\nNumber of Initial Conditions: {len(initial_conditions_pr)}\n")

    return initial_conditions
