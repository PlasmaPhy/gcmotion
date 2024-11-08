import numpy as np

from gcmotion.utils._logger_setup import logger


def calcWNU_grid(
    cwp,
    theta: np.array,
    psi: np.array,
    contour_Phi: bool,
    units: str,
):
    r"""Returns a single value or a grid of the calculated Hamiltonian.

    Only to be called internally, by ``energy_contour()``..

    Args:
        theta (np.array): The :math:`\theta` values.
        psi (np.array): The :math:`\psi` values.
        contour_Phi (bool): Whether or not to add the electric potential term
            :math:`e\Phi`.
        units (str): The energy units.
    """
    # Get all needed attributes first
    mu = cwp.muNU.magnitude
    Pzeta = cwp.Pzeta0NU.magnitude
    qfactor = cwp.qfactor
    bfield = cwp.bfield
    efield = cwp.efield

    b, _, g = bfield.bigNU(psi, theta)
    psip = qfactor.psipNU(psi)

    rho = (Pzeta + psip) / g
    W = (1 / 2) * rho**2 * b**2 + mu * b  # Without Φ

    # Add Φ if asked
    if contour_Phi:
        Phi = efield.PhiNU(psi, theta)
        W += Phi  # all normalized
        logger.debug("\tContour: Adding Φ term to the contour.")

    return W
