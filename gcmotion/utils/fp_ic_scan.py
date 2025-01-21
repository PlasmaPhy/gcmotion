import numpy as np

from gcmotion.entities.profile import Profile
from gcmotion.utils.fp_system import system  # Assuming this is available


def _find_local_minima(arr: np.ndarray):
    """
    Find local minima in a 2D array
    """

    nx, ny = arr.shape
    indices = []
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if (
                (arr[i, j] < arr[i + 1, j + 1])
                and (arr[i, j] < arr[i + 1, j])
                and (arr[i, j] < arr[i + 1, j - 1])
                and (arr[i, j] < arr[i - 1, j + 1])
                and (arr[i, j] < arr[i - 1, j])
                and (arr[i, j] < arr[i - 1, j - 1])
                and (arr[i, j] < arr[i, j + 1])
                and (arr[i, j] < arr[i, j - 1])
            ):

                indices.append((i, j))

    return indices


def fp_ic_scan(
    profile: Profile,
    method: str,
    theta_grid_density: int = 400,
    psi_grid_density: int = 400,
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.8],
    tol: float = 1e-7,
):
    theta_min, theta_max = theta_lim[0], 1.05 * theta_lim[1]
    psi_min, psi_max = psi_lim[0], 1.01 * psi_lim[1]

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
                system_values[i, j] = theta_dot**2 + (70 * psi_dot) ** 2

    # Find local minima
    indices = _find_local_minima(system_values)

    filtered_minima = [(i, j) for (i, j) in indices if system_values[i, j] < tol]

    # Extract minima positions
    fixed_point_candidates = [(theta_grid[i, j], psi_grid[i, j]) for i, j in filtered_minima]

    return fixed_point_candidates
