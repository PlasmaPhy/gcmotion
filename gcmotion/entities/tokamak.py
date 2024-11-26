from pint import Quantity
from termcolor import colored

from gcmotion.utils.logger_setup import logger
from gcmotion.tokamak.qfactor import QFactor
from gcmotion.tokamak.bfield import MagneticField
from gcmotion.tokamak.efield import ElectricField


class Tokamak:
    r"""Returns a Tokamak entity.

    A Tokamak entity holds all information about a device's dimensions, as well
    as the q-factor, and the magnetic and electric field inside it.

    Parameters
    ----------
    R : Quantity
        The tokamak's major radius dimensions of [length].
    a : Quantity
        The tokamak's minor radius dimensions of [length].
    qfactor : :py:class:`~gcmotion.qfactor.QFactor`
        Qfactor object that supports query methods for getting
        values of :math:`q(\psi)` and :math:`\psi_p(\psi)`.
    bfield : :py:class:`~gcmotion.bfield.MagneticField`
        Magnetic Field Object that supports query methods for
        getting values of the field magnitude, plasma currents and
        their derivatives.
    efield : :py:class:`~gcmotion.efield.ElectricField`
        Electric Field Object that supports query methods for
        getting values of the field itself and the derivatives of
        its potential.

    Notes
    -----
    A tokamak object contains the following attributes, which include the input
    arguements:

        #. R, RNU : Quantities
            The tokamak's major radius in [meters]/[NUmeters].
        #. a, aNU : Quantities
            The tokamak's minor radius in [meters]/[NUmeters].
        #. B0, B0NU : Quantities
            The strength of the magnetic field on the magnetic field in
            [Tesla]/[NUTesla].
        #. psi_wall, psi_wallNU : Quantities
            The :math:`\psi` value in the tokamak's wall.
        #. psip_wall, psip_wallNU : Quantities
            The :math:`\psi_p` value in the tokamak's wall.

    Example
    -------
    How to initialize a PhysicalParameters object:

    >>> import gcmotion as gcm
    >>>
    >>> # Quantity Constructor
    >>> Rnum = 1.65
    >>> anum = 0.5
    >>> B0num = 1
    >>> species = "p"
    >>> ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)
    >>>
    >>> # Intermediate values
    >>> R = Q(Rnum, "meters")
    >>> a = Q(anum, "meters")
    >>> B0 = Q(B0num, "Tesla")
    >>> i = Q(0, "NUPlasma_current")
    >>> g = Q(1, "NUPlasma_current")
    >>> Ea = Q(73500, "Volts/meter")
    >>>
    >>> tokamak = gcm.Tokamak(
    ...     R=R,
    ...     a=a,
    ...     qfactor=gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
    ...     bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
    ...     efield=gcm.efield.Radial(a, Ea, B0, peak=0.98, rw=1 / 50),
    ... )

    """

    def __init__(
        self,
        R: Quantity,
        a: Quantity,
        qfactor: QFactor,
        bfield: MagneticField,
        efield: ElectricField,
    ):
        r"""Sets up fields in the correct units."""
        logger.info("==> Initializing Tokamak")

        # Construct objects
        self.qfactor = qfactor
        self.bfield = bfield
        self.efield = efield

        # Construct Quantities
        self.R = R.to("meters")
        self.a = a.to("meters")
        self.B0 = bfield.B0.to("Tesla")
        self.RNU = self.R.to("NUmeters")
        self.aNU = self.a.to("NUmeters")
        self.B0NU = self.B0.to("NUTesla")

        # Calculate last closed surfaces
        self.psi_wall = (self.B0 * self.a**2 / 2).to("Magnetic_flux")
        self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")
        self.psip_wallNU = (
            self.qfactor.psipNU(self.psi_wallNU.magnitude)
            * self.psi_wallNU.units
        )
        self.psip_wall = self.psip_wallNU.to("Magnetic_flux")

        logger.info("\tTokamak Initialization Complete")

    def __repr__(self):
        return (
            "Tokamak:"
            + f"R = {self.R:.4g~}, "
            + f"a = {self.a:.4g~}, "
            + f"qfactor = [{self.qfactor}], "
            + f"bfield = [{self.bfield}], "
            + f"efield = [{self.efield}]"
        )

    def __str__(self):
        return (
            colored("\nTokamak:\n", "green")
            + f"{'R':>23} : {f'{self.R:.4g~}':<16}({self.RNU:.4g~})\n"
            + f"{'a':>23} : {f'{self.a:.4g~}':<16}({self.aNU:.4g~})\n"
            + f"{'B0':>23} : {f'{self.B0:.4g~}':<16}({self.B0NU:.4g~})\n"
            + f"{'q-factor':>23} : {self.qfactor}\n"
            + f"{'Magnetic Field':>23} : {self.bfield}\n"
            + f"{'Electric Field':>23} : {self.efield}\n"
            + f"{'psi_wall':>23} : {f'{self.psi_wall:.4g~}':<16}"
            + f"({self.psi_wallNU:.4g~})\n"
            + f"{'psip_wall':>23} : {f'{self.psip_wall:.4g~}':<16}"
            + f"({self.psip_wallNU:.4g~})"
        )
