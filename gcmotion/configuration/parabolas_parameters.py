from dataclasses import dataclass
import numpy as np

figsize = 13, 7  # Global window size
dpi = 100  # Global dpi


@dataclass()
class ParabolasConfig:
    Pzetalim: tuple = (-1.5, 1)  # result after division by psip_wall.m
    Pzeta_density: int = 1000


@dataclass
class ParabolasPlotConfig:
    # Figure keywords
    figsize: tuple = figsize
    dpi: int = dpi
    title: str = (
        r"Constant-$\mu$ slices (plane cuts) of the three dimensional COM space (E,$\mu,P_{\zeta}$)"
    )
    layout: str = "constrained"
    facecolor: str = "lightskyblue"
    linewidth: int = 2
    xlabel_fontsize: int = 13
    ylabel_fontsize: int = 13
    ylabel_rotation: int = 0
    legend: bool = True
    # Parabolas keywords
    enlim: tuple = (0, 3)
    Pzetalim: tuple = (-1, 1)  # result after division by psip_wall.m
    Pzeta_density: int = 1000
    LAR_TPB: bool = False
