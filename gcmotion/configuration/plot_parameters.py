from dataclasses import dataclass
from numpy import pi

figsize = 16, 9  # Global window size
dpi = 90  # Global dpi


@dataclass()
class ProfileEnergyContourConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "lightskyblue"
    # Default parameter values
    thetalim: tuple = (-pi, pi)
    psilim: tuple = (0, 1.2)  # times psi_wall
    levels: int = 30
    E_units: str = "keV"
    flux_units: str = "Tesla * meter^2"
    potential: bool = True
    wall: bool = True
    cursor: bool = True  # Mild performance hit
    # Colorbar
    numticks: int = 10
    cbarlabelsize: int = 12


@dataclass()
class ProfilePzetaContourConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "lightskyblue"
    # Default parameter values
    zetalim: tuple = (-pi, pi)
    psilim: tuple = (0, 1.2)  # times psi_wall
    levels: int = 30
    Pzeta_units: str = "Tesla * meter^2"
    flux_units: str = "Tesla * meter^2"
    potential: bool = True
    wall: bool = True
    cursor: bool = True  # Mild performance hit
    # Colorbar
    numticks: int = 10
    cbarlabelsize: int = 12


@dataclass()
class ParticleEvolutionConfig:
    # Default parameter values
    which: str = "all"
    units: str = "SI"
    percentage: int = 100
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "lightskyblue"
    titlesize: float = 20
    titlecolor: str = "blue"
    # Default parameter values
    which: str = "all"
    units: str = "SI"
    percentage: int = 100
    # Scatter kw
    s: float = 0.2
    color: str = "blue"
    marker: str = "o"
    labelsize: int = 10
    labelpad: float = 8


@dataclass
class QfactorProfileConfig:
    # Figure keywords
    figsize: tuple = (16, 5)
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "lightskyblue"
    titlesize: float = 20
    titlecolor: str = "blue"
    # Default parameter values
    span: tuple = (0, 1.1)
    norm: bool = True
    show: bool = True
    # Plot options
    points: int = 200
    wall_color: str = "red"
    qwall_color: str = "black"
    qwall_style: str = "--"
    psip_wall_color: str = "black"
    psip_wall_style: str = "--"
    labelsize: float = 10
    ax_title_size: float = 20
