import numpy as np

from gcmotion.entities.profile import Profile


def energy_Ptheta(
    psi: float | np.ndarray,
    theta: float | np.ndarray,
    profile: Profile,
    contour_Phi: bool,
):
    r"""Calculates the Energy and :math:`P_\theta` as a function of
    :math:`\psi, \theta, \mu, P_\zeta` and the tokamak's
    profile in [NU]..

    Only to be called internally, by
    :py:func:`~gcmotion.plotters.energy_contour.energy_contour()`..

    Parameters
    ----------
    psi : float | np.array
        The :math:`\psi` values in [NU]
    theta : float | np.array
        The :math:`\theta` values.
    profile : Profile
        The tokamak Rrofile.
    contour_Phi : bool
        Whether or not to add the electric potential term :math:`\Phi`.

    Returns
    -------
    tuple of (float | np.ndarray) (depending in the input)
        The calculated Energy values in [NU] and Ptheta in [NU].

    """
    # Get all needed attributes first
    mu = profile.muNU.m
    Pzeta = profile.PzetaNU.m
    qfactor = profile.qfactor
    bfield = profile.bfield
    efield = profile.efield

    b, i, g = bfield.bigNU(psi, theta)

    # Calculate psi and psip, in NU
    psip = qfactor.psipNU(psi)

    rho = (Pzeta + psip) / g
    W = (1 / 2) * rho**2 * b**2 + mu * b  # Without Φ

    # Add Φ if asked
    if contour_Phi:
        Phi = efield.PhiNU(psi, theta)
        W += Phi  # all normalized

    Ptheta = psi + rho * i

    return W, Ptheta
