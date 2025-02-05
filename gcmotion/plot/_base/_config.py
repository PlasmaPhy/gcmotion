r"""
Sets up configuration for base plotting scripts.

All these values can be overwritten if passed as an arguement to the
corresponding function.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class _ProfileEnergyContourConfig:
    # Default optional arguements
    thetalim: tuple = (-np.pi, np.pi)
    psilim: tuple = (0, 1.2)
    levels: int = 30
    flux_units: str = "Tesla * meters^2"
    E_units: str = "keV"
    potential: bool = True
    wall: bool = True
    # Contour
    projection: str | None = None  # None = default
    mode: str = "filled"  # "filled" or "lines"
    grid_density: int = 200
    cmap: str = "plasma"
    locator: str = "log"
    log_base: float = 1.0001
    zorder: int = 0
    Pthetaax: bool = True
    cursor: bool = True
    # Labels
    labelsize: float = 15
    ticknum: int = 10
    ticksize: float = 12


@dataclass
class _ProfilePzetaContourConfig:
    # Default optional arguements
    thetalim: tuple = (-np.pi, np.pi)
    psilim: tuple = (0, 1.2)
    levels: int = 30
    flux_units: str = "Tesla * meters^2"
    Pzeta_units: str = "Tesla * meters^2"
    potential: bool = True
    wall: bool = True
    # Contour
    mode: str = "filled"  # "filled" or "lines"
    grid_density: int = 200
    cmap: str = "plasma"
    locator: str = ""  # Pzeta can be negative so dont use LogLocator
    log_base: float = 1.0001
    zorder: int = 0
    Pthetaax: bool = True
    cursor: bool = True
    # Labels
    labelsize: float = 15
    ticknum: int = 10
    ticksize: float = 12


@dataclass
class _ColorbarConfig:
    location: str = "top"
    numticks: int = 15
    label: str = ""
    labelsize: float = 15


@dataclass
class _ParticlePoloidalDrift:
    units: str = "SI"  # "SI" or "NU"
    thetalim: tuple = (-np.pi, np.pi)
    percentage: int = 100
    # Scatter kw
    s: float = 1
    color: str = "r"
    marker: str = "."
