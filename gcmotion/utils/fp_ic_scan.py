import numpy as np
import pint
from collections import namedtuple, deque
from functools import lru_cache
from numba import njit

from gcmotion.utils.fp_system import system  # Assuming this is available

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


# Cached version of system
@lru_cache(maxsize=None)
def cached_system(theta, psi, parameters, profile):
    return system(theta, psi, parameters=parameters, profile=profile)


@njit
def evaluate_neighbors(system_values, i, j):
    """Helper function to evaluate neighbors of a grid point."""
    neighbors = [
        system_values[i - 1, j - 1],
        system_values[i - 1, j],
        system_values[i - 1, j + 1],
        system_values[i, j - 1],
        system_values[i, j + 1],
        system_values[i + 1, j - 1],
        system_values[i + 1, j],
        system_values[i + 1, j + 1],
    ]
    return neighbors


def fp_ic_scan(
    parameters: namedtuple,
    profile: namedtuple,
    theta_grid_density: int = 800,
    psi_grid_density: int = 800,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.8],
    tol: float = 1e-6,
):
    # Define theta and psi bounds
    theta_min, theta_max = theta_lim
    psi_min, psi_max = 0.99 * psi_lim[0], 1.01 * psi_lim[1]

    # Create the grid
    theta_grid, psi_grid = np.meshgrid(
        np.linspace(theta_min, theta_max, theta_grid_density),
        np.linspace(psi_min, psi_max, psi_grid_density),
    )

    grid_shape = theta_grid.shape

    # Precompute all system values for the grid
    system_values = np.empty_like(theta_grid)
    for i in range(grid_shape[0]):
        for j in range(grid_shape[1]):
            system_values[i, j] = cached_system(
                theta_grid[i, j], psi_grid[i, j], parameters, profile
            )

    # Find fixed_point_candidates (minima/maxima/saddle points) using a combination
    # of neighbor checks
    fixed_point_candidates = deque([])
    for i in range(1, grid_shape[0] - 1):  # Exclude first and last row
        for j in range(1, grid_shape[1] - 1):  # Exclude first and last column
            center = system_values[i, j]

            # Check if center value is below tolerance
            if center >= tol:
                continue

            # Check neighbors
            neighbors = evaluate_neighbors(system_values, i, j)

            # If center is smaller than all neighbors, add to fixed_point_candidates
            if all(center < n for n in neighbors):
                fixed_point_candidates.append((theta_grid[i, j], psi_grid[i, j]))

    # Convert fixed_point_candidates to list of initial conditions
    initial_conditions = list(fixed_point_candidates)

    return initial_conditions