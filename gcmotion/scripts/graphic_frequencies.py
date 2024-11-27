r"""
Calculates the Energy-frequency relations from the areas under the contour
plots.
"""

import numpy as np
import matplotlib.pyplot as plt

from gcmotion.utils.energy_Ptheta import energy_Ptheta


def graphic_frequencies(psi_lim: list, parameters: dict):

    # Unpack Parameters
    mu = parameters["NU"]
    Pzeta = parameters["Pzeta"]
    tokamak = parameters["tokamak"]
    Q = parameters["Q"]

    # Set Quantities
    R = tokamak["R"]
    a = tokamak["a"]
    B0 = tokamak["B0"]
    qfactor = tokamak["qfactor"]
    bfield = tokamak["bfield"]
    efield = tokamak["efield"]

    theta_span = [0, 2 * np.pi]
    psi_span = Q(psi_lim, "NUpsi_wall")
    print(psi_lim)

    theta_grid, psi_grid = np.meshgrid(
        np.linspace(theta_span[0], theta_span[1], 100),
        np.linspace(psi_span[0], psi_span[1], 100),
    )

    E_grid, Ptheta_grid = energy_Ptheta(
        psi_grid, theta_grid, mu, Pzeta, profile=profile, contour_Phi=True
    )
