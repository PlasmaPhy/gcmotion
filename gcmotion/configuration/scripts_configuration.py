import numpy as np
from dataclasses import dataclass


@dataclass
class SolverConfig:
    atol: float = 1e-12  # Scipy's default is 1e-6
    rtol: float = 1e-12  # Scipy's default is 1e-3


@dataclass
class FrequencyConfig:
    # Figure parameters
    figsize: tuple = (16, 9)
    dpi: int = 90
    layout: str = "constrained"
    # Resulting mosaic
    mosaic: str = "debug"

    # Auto
    auto_thetalim: tuple = (-np.pi, np.pi)
    auto_psilim: tuple = (0, 1.2)
    auto_levels: int = 100

    # Energy contour (omega_theta)
    psilim: tuple = (0, 1.2)
    levels: int = 250
    grid_density: int = 200
    # Misc contour plot params
    flux_units: str = "Tesla * meters^2"
    E_units: str = "keV"
    potential: bool = True
    wall: bool = True
    st_end_points: bool = False
    # Filters
    theta_rtol: float = (4 * np.pi) / 100
    psi_rtol: float = 1 / 100


@dataclass
class QKineticConfig:
    levels: int = 100
