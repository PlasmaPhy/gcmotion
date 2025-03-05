from dataclasses import dataclass
from numpy import pi
import numpy as np

figsize = 13, 7  # Global window size
dpi = 100  # Global dpi


@dataclass()
class ProfileEnergyContourConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    projection: str | None = None  # None = default
    # Default parameter values
    thetalim: tuple = (-pi, pi)
    psilim: tuple = (0, 1.2)  # times psi_wall
    levels: int = 30
    E_units: str = "keV"
    flux_units: str = "Tesla * meter^2"
    canmon_units: str = "NUCanonical_momentum"
    potential: bool = True
    wall: bool = True
    cursor: bool = True  # Mild performance hit
    show: bool = True
    # Colorbar
    numticks: int = 10
    cbarlabelsize: int = 12


@dataclass()
class ProfilePzetaContourConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    # Default parameter values
    zetalim: tuple = (-pi, pi)
    psilim: tuple = (0, 1.2)  # times psi_wall
    levels: int = 30
    Pzeta_units: str = "Tesla * meter^2"
    flux_units: str = "Tesla * meter^2"
    potential: bool = True
    wall: bool = True
    cursor: bool = True  # Mild performance hit
    show: bool = True
    # Colorbar
    numticks: int = 10
    cbarlabelsize: int = 12


@dataclass
class QfactorProfileConfig:
    # Figure keywords
    figsize: tuple = (13, 5)
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    titlesize: float = 20
    titlecolor: str = "blue"
    # Default parameter values
    span: tuple = (0, 1.1)
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


@dataclass
class EfieldProfileConfig:
    # Figure keywords
    figsize: tuple = (13, 5)
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    titlesize: float = 20
    titlecolor: str = "blue"
    # Default parameter values
    span: tuple = (0, 1.1)
    show: bool = True
    # Plot options
    points: int = 400
    wall_color: str = "red"
    labelsize: float = 10
    ax_title_size: float = 20
    # Units options
    field_units: str = "kV/m"
    potential_units: str = "kV"


@dataclass
class MagneticProfileConfig:
    # Figure keywords
    figsize: tuple = (13, 7)
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    titlesize: float = 20
    titlecolor: str = "blue"
    # Default parameter values
    span: tuple = (0, 1.1)
    units: str = "NU"
    plot_derivatives: bool = True
    coord: str = "r"  # psi / r
    show: bool = True
    # Contour options
    grid_density: int = 200
    levels: int = 20
    bcmap: str = "inferno"
    icmap: str = "viridis"
    gcmap: str = "viridis"
    locator: str = ""
    log_base: float = 1.00001
    # 2d plot options
    current_color: str = "b"
    derivative_color: str = "r"
    linewidth: float = 2
    # Label options
    labelsize: float = 15
    ax_title_size: float = 20
    ax_title_pad: float = 25


@dataclass()
class ParticleEvolutionConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    titlesize: float = 20
    titlecolor: str = "blue"
    # Default parameter values
    which: str = "all"
    units: str = "SI"
    percentage: int = 100
    show: bool = True
    # Scatter kw
    s: float = 0.2
    color: str = "blue"
    marker: str = "o"
    labelsize: int = 10
    labelpad: float = 8


@dataclass()
class ParticlePoloidalDrift:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    # Default parameter values
    projection: str | None = None  # None = default
    thetalim: tuple = (-pi, pi)
    psilim: str | tuple = "auto"  # times psi_wall, or "auto"
    levels: int = 30
    E_units: str = "keV"
    flux_units: str = "Tesla * meter^2"
    potential: bool = True
    initial: bool = True
    wall: bool = True
    cursor: bool = True  # Mild performance hit
    show: bool = True
    # Colorbar
    numticks: int = 10
    cbarlabelsize: int = 12


@dataclass
class AutoYspan:
    zoomout: float = 0.75
    hardylim: float = 3  # times psi_wall


@dataclass
class ParabolasPlotConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    linewidth: int = 2
    # Title keywords
    title_fontsize: float = 15
    title_color: str = "black"
    # Labels keywords
    xlabel_fontsize: float = 13
    xlabel_rotation: int = 0
    ylabel_fontsize: float = 13
    ylabel_rotation: int = 0
    # Legend keywords
    parabolas_legend: bool = True
    # Parabolas keywords
    enlim: tuple = (0, 3)
    Pzetalim: tuple = (-1, 1)  # result after division by psip_wall.m
    Pzeta_density: int = 1000
    TPB_density: int = 100
    plot_TPB: bool = False
    parabolas_color: str = "orange"
    TPB_X_color: str = "#E65100"
    TPB_O_color: str = "#1f77b4"
    TPB_X_linestyle: str = "solid"
    TPB_O_linestyle: str = "solid"
    TPB_X_markersize: float = 2
    TPB_O_markersize: float = 2
    # Dashed line keywords
    show_d_line: bool = True
    d_line_color: str = "black"
    d_linewidth: int = 1
    d_line_alplha: float = 0.5


@dataclass()
class BifurcationPlotConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    sharex: bool = True
    # x Label keywords
    xlabel_fontsize: float = 10
    # Suptitle keywords
    suptitle_fontsize: float = 15
    suptitle_color: str = "black"
    # Bifurcation keywords
    thetalim: tuple = (-np.pi, np.pi)
    psilim: tuple = (0, 1.8)
    plot_energy_bif: bool = True
    which_COM: str = "Pzeta"


@dataclass
class RZContoursConfig:
    # Figure keywords
    figsize: tuple = (6, 8)
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    show: bool = True
    # Contour keywords
    cmap: str = "plasma"
    levels: int = 20
    mode: str = None
    units: str = "NUmf"
    which_Q: str = "flux"
    # Locator keywords
    log_base: float = 1.0001
    locator: str = "linear"
    # Boundary keywords
    black_boundary: bool = True
    boundary_linewidth: int = 1
    # Stationary curves keywords
    plot_stationary_curves: bool = True
    stat_curves_color: str = "black"
    stat_curves_linewidth: float = 1
    stat_curves_linestyle: str = "dashed"
    # Labels - Title keywords
    xlabel_fontsize: float = 15
    ylabel_fontsize: float = 15
    title_fontsize: float = 15
    # Colorbar keywords
    cbarlabel_fontsize: float = 10
    cbar_ticks: int = 10
    # Numerical keywords
    parametric_density: int = 500
    xmargin_perc: float = 0.1
    ymargin_perc: float = 0.1


@dataclass
class RZBigContoursConfig:
    # Figure keywords
    figsize_B: tuple = (17, 8)
    figsize_I: tuple = (11, 8)
    figsize_g: tuple = (11, 8)
    dpi: int = dpi
    layout: str = "constrained"
    facecolor: str = "white"
    # B figure keywords
    B_suptitle_fontsize: float = 15
    B_suptitle_color: str = "black"
    B_units: str = "Tesla"
    # I figure keywords
    I_suptitle_fontsize: float = 15
    I_suptitle_color: str = "black"
    I_units: str = "NUpc"
    # g figure keywords
    g_suptitle_fontsize: float = 15
    g_suptitle_color: str = "black"
    g_units: str = "NUpc"
    # Numerical keywords
    parametric_density: int = 500
