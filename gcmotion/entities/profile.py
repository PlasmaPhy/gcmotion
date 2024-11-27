r"""
==============
Profile Entity
==============

This module defines the "Profile" entity, which is a child class of the
"Tokamak" and "PhysicalParameters" classes.

A lot of analysis can be done upon this class, since it essentially specifies a
family of particles with the same 3 Constants of Motion, E, mu and Pzeta, in a
specific tokamak device, which fully define their Hamiltonian.

This class also constructs the QuantityConstructor to be used internally, so
every subclass should grab it from here instead of redifining it.
"""

from gcmotion.utils.logger_setup import logger
from gcmotion.utils.quantity_constructor import QuantityConstructor

from gcmotion.entities.tokamak import Tokamak
from gcmotion.entities.physical_parameters import PhysicalParameters


class Profile(Tokamak, PhysicalParameters):
    r"""Returns a Profile entity.

    A Profile entity represents the subset of all Particles of certain species,
    :math:`\mu` and :math:`P_\zeta` in a specific tokamak configuration.

    Parameters
    ----------
    tokamak : Tokamak
        The Tokamak entity.
    params : PhysicalParameters
        The PhysicalParameters entity.

    Example
    -------
    How to create a `Profile` object.

    >>> import gcmotion as gcm
    >>>
    >>> #Quantity Constructor
    >>> Rnum = 1.65
    >>> anum = 0.5
    >>> B0num = 1
    >>> species = "p"
    >>> Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)
    >>>
    >>> # Intermediate Quantities
    >>> R = Q(Rnum, "meters")
    >>> a = Q(anum, "meters")
    >>> B0 = Q(B0num, "Tesla")
    >>> i = Q(0, "NUPlasma_current")
    >>> g = Q(1, "NUPlasma_current")
    >>> Ea = Q(73500, "Volts/meter")
    >>>
    >>> # Construct a Tokamak
    >>> tokamak = gcm.Tokamak(
    ...     R=R,
    ...     a=a,
    ...     qfactor=gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2),
    ...     bfield=gcm.bfield.LAR(B0=B0, i=i, g=g),
    ...     efield=gcm.efield.Radial(a, Ea, B0, peak=0.98, rw=1 / 50),
    ... )
    >>>
    >>> # Construct a PhysicalParameters set
    >>> params = gcm.PhysicalParameters(
    ...     species=species,
    ...     mu=Q(1e-5, "NUMagnetic_moment"),
    ...     Pzeta=Q(-0.0272, "NUMagnetic_flux"),
    ... )
    >>>
    >>> # Create a Profile
    >>> profile = gcm.Profile(tokamak, params)

    .. admonition:: For Developers

        The Profile now has all the information to initialize the Quantity
        Constructor (again) to be used internally. All child classes should
        grab it from here.

    """

    def __init__(self, tokamak: Tokamak, params: PhysicalParameters):

        logger.info("==> Initializing Profile...")

        # Grab attributes from parents
        Tokamak.__init__(
            self,
            R=tokamak.R,
            a=tokamak.a,
            qfactor=tokamak.qfactor,
            bfield=tokamak.bfield,
            efield=tokamak.efield,
        )

        PhysicalParameters.__init__(
            self,
            species=params.species,
            mu=params.mu,
            Pzeta=params.Pzeta,
        )

        # Store those for easier reference
        self.tokamak = tokamak
        self.params = params

        # This is a good place to initialize Q, since we have everything we
        # need. Every child class should grab it from here.
        self.Q = QuantityConstructor(
            R=self.R.magnitude,
            a=self.a.magnitude,
            B0=self.B0.magnitude,
            species=self.species,
        )

        logger.info("\tProfile setup complete.")

    def __repr__(self):
        return PhysicalParameters.__repr__(self) + Tokamak.__repr__(self)

    def __str__(self):
        return PhysicalParameters.__str__(self) + Tokamak.__str__(self)
