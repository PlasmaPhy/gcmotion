import numpy as np

from gcmotion.entities.profile import Profile
from gcmotion.utils.fp_system import system  # Assuming this is available


def fp_ic_scan(
    profile: Profile,
    method: str,
    theta_grid_density: int = 400,
    psi_grid_density: int = 400,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.8],
    tol: float = 1e-7,
):
    theta_min, theta_max = theta_lim
    psi_min, psi_max = 0.99 * psi_lim[0], 1.01 * psi_lim[1]

    theta_grid, psi_grid = np.meshgrid(
        np.linspace(theta_min, theta_max, theta_grid_density),
        np.linspace(psi_min, psi_max, psi_grid_density),
    )

    system_values = np.empty_like(theta_grid)
    grid_shape = theta_grid.shape

    # Precompute system values
    for i in range(grid_shape[0]):
        for j in range(grid_shape[1]):
            if method == "differential evolution":
                system_values[i, j] = system(
                    theta_grid[i, j], psi_grid[i, j], profile, method=method
                )
            elif method == "fsolve":
                theta_dot, psi_dot = system(
                    theta_grid[i, j], psi_grid[i, j], profile, method=method
                )
                system_values[i, j] = theta_dot**2 + psi_dot**2

    # Find local minima
    padded_values = np.pad(system_values, pad_width=1, mode="edge")
    neighbors = np.stack(
        [
            padded_values[1:-1, :-2],
            padded_values[1:-1, 2:],  # Left, right
            padded_values[:-2, 1:-1],
            padded_values[2:, 1:-1],  # Up, down
            padded_values[:-2, :-2],
            padded_values[:-2, 2:],  # Diagonals
            padded_values[2:, :-2],
            padded_values[2:, 2:],
        ],
        axis=0,
    )

    is_minima = (system_values < np.min(neighbors, axis=0)) & (system_values < tol)

    # Extract minima positions
    minima_indices = np.argwhere(is_minima)
    fixed_point_candidates = [(theta_grid[i, j], psi_grid[i, j]) for i, j in minima_indices]

    return fixed_point_candidates
