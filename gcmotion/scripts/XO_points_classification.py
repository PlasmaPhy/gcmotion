import numpy as np
from collections import namedtuple, deque
from gcmotion.utils.energy_Ptheta import energy_Ptheta


def XO_points_classification(
    unclassified_fixed_points: np.ndarray,
    parameters: namedtuple,
    profile: namedtuple,
    grid_density: int = 10,
):
    # Parameters
    mu = parameters.mu
    Pzeta0 = parameters.Pzeta0

    O_points = deque()  # Deque for stable O-points
    X_points = deque()  # Deque for unstable X-points and saddle X-points

    for fixed_point in unclassified_fixed_points:

        theta_fixed = fixed_point[0]
        P_theta_fixed = fixed_point[1]

        # Define the grid around the fixed point
        theta_min = 0.9 * theta_fixed
        theta_max = 1.1 * theta_fixed

        if theta_fixed < 0:
            theta_min = 1.1 * theta_fixed
            theta_max = 0.9 * theta_fixed
        elif theta_fixed < 1e-2:
            theta_min = theta_fixed - 0.1
            theta_max = theta_fixed + 0.1

        P_theta_min = 0.9 * P_theta_fixed
        P_theta_max = 1.1 * P_theta_fixed

        # Create the grid
        theta, P_thetaNU = np.meshgrid(
            np.linspace(theta_min, theta_max, grid_density),
            np.linspace(P_theta_min, P_theta_max, grid_density),
        )

        # Compute the Hamiltonian values on the grid
        WNU, _ = energy_Ptheta(
            psi=P_thetaNU, theta=theta, Pzeta=Pzeta0, mu=mu, profile=profile, contour_Phi=True
        )  # BECAUSE psi=P_thetaNU ONLY WORKS FOR LAR AT THE MOMENT

        # Locate the grid index of the fixed point
        idx_theta = np.abs(theta[0, :] - theta_fixed).argmin()  # Closest row in theta
        idx_P_theta = np.abs(P_thetaNU[:, 0] - P_theta_fixed).argmin()  # Closest column in psi

        # Energy at the fixed point
        W_fixed = WNU[idx_theta, idx_P_theta]

        # Extract the energy of the 4 second closest neighbors to avoid
        # numerical errors
        neighbor_up = WNU[idx_theta, idx_P_theta + 2]
        neighbor_down = WNU[idx_theta, idx_P_theta - 2]
        neighbor_right = WNU[idx_theta + 2, idx_P_theta]
        neighbor_left = WNU[idx_theta - 2, idx_P_theta]

        # Classify the fixed point
        if all(W > W_fixed for W in [neighbor_up, neighbor_down, neighbor_left, neighbor_right]):
            O_points.append([theta_fixed, P_theta_fixed])  # Add to O-points deque
        elif all(W < W_fixed for W in [neighbor_up, neighbor_down, neighbor_left, neighbor_right]):
            X_points.append([theta_fixed, P_theta_fixed])  # Add to X-points deque (unstable max)
        elif (
            neighbor_up < W_fixed
            and neighbor_down < W_fixed
            and neighbor_left > W_fixed
            and neighbor_right > W_fixed
        ) or (
            neighbor_up > W_fixed
            and neighbor_down > W_fixed
            and neighbor_left < W_fixed
            and neighbor_right < W_fixed
        ):
            X_points.append([theta_fixed, P_theta_fixed])  # Add to X-points deque (saddle)

    return X_points, O_points
