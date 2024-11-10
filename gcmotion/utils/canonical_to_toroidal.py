import numpy as np

from gcmotion.utils.logger_setup import logger


def canonical_to_toroidal(
    cwp, percentage: int = 100, truescale: bool = True
) -> tuple:
    r"""Calculates the toroidal coordionates of the particles orbit,
    :math:`(r, \theta, \zeta)`.

    :math:`r = \sqrt{2\psi}` rather than :math:`\psi` itself is used for
    the plot, since it is a better representation of the actual orbit.

    :meta public:

    Parameters
    ----------
    percentage : int, optional
        The percentage of the orbit to be plotted. Defaults to 100.

    truescale : bool, optional
        Whether or not to use the actual tokamak dimensions, or fit
        them around the orbit for better visibility. Defaults to True.

    Returns
    -------
    5-tuple of 2 floats and 3 np.arrays:
        The major and minor radii of the (possibly scaled) tokamak and
        the toroidal coordionates of the particles orbit
        :math:`(r, \theta, \zeta)`.

    """
    logger.info("Calculating torus plotting points...")

    # Get all needed attributes first
    R, a = cwp.R.magnitude, cwp.a.magnitude
    theta = cwp.theta.magnitude
    psi = cwp.psi.magnitude
    zeta = cwp.zeta.magnitude

    if percentage < 1 or percentage > 100:
        percentage = 100
        print("Invalid percentage. Plotting the whole thing.")
        logger.warning("Invalid percentage: Plotting the whole thing...")

    points = int(np.floor(theta.shape[0] * percentage / 100) - 1)
    theta_torus = theta[:points]
    z_torus = zeta[:points]
    r_torus = np.sqrt(2 * psi[:points])

    # Torus shape parameters
    r_span = [r_torus.min(), r_torus.max()]
    logger.debug(
        f"\tr-span calculated:[{r_span[0]:.4g}, {r_span[1]:.4g}]m, with a={a}m."
    )

    if truescale:
        Rtorus = R
        atorus = a
        logger.debug("\tReturning toroidal coordinates in True scale.")
    else:
        Rtorus = (r_span[1] + r_span[0]) / 2
        atorus = 1.1 * Rtorus / 2
        r_torus *= 1 / 2
        logger.warning("Returning toroidal coordinates 'zoomed' in.")

    logger.info("--> Torus points calculation successful.")

    return Rtorus, atorus, r_torus, theta_torus, z_torus
