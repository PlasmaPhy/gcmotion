from dataclasses import dataclass


@dataclass(frozen=True, repr=True)
class SolverConfig:
    atol: float = 1e-8  # Scipy's default is 1e-6
    rtol: float = 1e-8  # Scipy's default is 1e-3


@dataclass(frozen=True, repr=True)
class FrequencyConfig:
    # Figure parameters
    figsize: tuple = (16, 9)
    dpi: int = 90
    layout: str = "constrained"
    # ψ limit and levels for countour lines search
    psilim: tuple = (0, 1.2)
    levels: int = 50
    grid_density: int = 200
    # Misc contour plot params
    flux_units: str = "Tesla * meters^2"
    E_units: str = "keV"
    potential: bool = True
    wall: bool = True
