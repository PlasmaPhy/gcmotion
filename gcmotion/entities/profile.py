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

import pint

from gcmotion.utils.logger_setup import logger
from gcmotion.utils.quantity_constructor import QuantityConstructor

from gcmotion.entities.tokamak import Tokamak
from gcmotion.entities.physical_parameters import PhysicalParameters
from gcmotion.configuration.physical_constants import PhysicalConstants

type Quantity = pint.Quantity


class Profile(Tokamak, PhysicalParameters):
    r"""

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
            E=params.E,
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

        # Calculate mass and charge
        _miNU = getattr(PhysicalConstants, self.species + "_M", None)
        _qiNU = getattr(PhysicalConstants, self.species + "_Z", None)
        self.miNU = self.Q(_miNU, "Proton_mass")
        self.qiNU = self.Q(_qiNU, "Proton_charge")
        self.mi = self.miNU.to("kilogram")
        self.qi = self.qiNU.to("Coulomb")

        logger.info("\tProfile setup complete.")

    def findPtheta(self, psi: Quantity):
        r"""Calculates Ptheta from psi. Output units are the same as input
        units.

        Only applicable in the absence of perturbations.

        Parameters
        ----------
        psi : Quantity
            The psi Quantity.

        Returns
        -------
        Quantity
            The calculated Ptheta Quantity.
        """
        # Store input units
        input_units = psi.units

        # All Quantities and operations are in NU, and the conversion to input
        # units takes place at the end.
        psiNU = psi.to("NUMagnetic_flux")
        _psiNU = psiNU.magnitude

        # Calculate and Quantify currents
        _, _iNU, _gNU = self.bfield.bigNU(_psiNU, 0)
        iNU = self.Q(_iNU, "NUPlasma_current")
        gNU = self.Q(_gNU, "NUPlasma_current")

        # Calculate and Quantify psip
        _psipNU = self.qfactor.psipNU(_psiNU)
        psipNU = self.Q(_psipNU, "NUmagnetic_flux")

        rhoNU = (self.PzetaNU + psipNU) / gNU
        PthetaNU = psiNU + rhoNU * iNU

        # Convert to input units and return
        return PthetaNU.to(input_units)

    def findEnergy(
        self, psi: Quantity, theta: float, units: str, potential: bool = True
    ):
        r"""Calculates the Energy of a particle characterized by a (psi, theta)
        pair.

        Parameters
        ----------
        psi : Quantity
            The particle's psi Quantity.
        theta : float
            The particle's :math:`\theta` angle.
        units : str
            The returned Energy units.
        potential : bool
            Whether or not to add the electric potential term in the energy.
            Defaults to True.

        Returns
        -------
        Quantity
            The calculated Energy Quantity.
        """
        # All Quantities and operations are in NU, and the conversion to input
        # units takes place at the end.
        psiNU = psi.to("NUMagnetic_flux")
        _psiNU = psiNU.magnitude

        # Calculate and Quantify bfield and currents
        _bNU, _iNU, _gNU = self.bfield.bigNU(_psiNU, theta)
        bNU = self.Q(_bNU, "NUTesla")
        gNU = self.Q(_gNU, "NUPlasma_current")

        # Calculate and Quantify psip
        _psipNU = self.qfactor.psipNU(_psiNU)
        psipNU = self.Q(_psipNU, "NUmagnetic_flux")

        rhoNU = (self.PzetaNU + psipNU) / gNU

        EnergyNU = (  # WARN: Unsure if q and m must appear here
            self.qiNU**2 / (2 * self.miNU)
        ) * rhoNU**2 * bNU**2 + self.muNU * bNU  # Without potential

        if potential:
            _PhiNU = self.efield.PhiNU(_psiNU, theta)
            PhiNU = self.Q(_PhiNU, "NUVolts")
            EnergyNU += self.qiNU * PhiNU

        # Convert to input units and return
        return EnergyNU.to(units)

    def findPzeta(
        self, psi: Quantity, theta: float, units: str, potential: bool = True
    ):

        # All Quantities and operations are in NU, and the conversion to input
        # units takes place at the end.
        psiNU = psi.to("NUMagnetic_flux")
        _psiNU = psiNU.magnitude

        # Calculate and Quantify bfield and currents
        _bNU, _iNU, _gNU = self.bfield.bigNU(_psiNU, theta)
        bNU = self.Q(_bNU, "NUTesla")
        gNU = self.Q(_gNU, "NUPlasma_current")

        # Calculate and Quantify psip
        _psipNU = self.qfactor.psipNU(_psiNU)
        psipNU = self.Q(_psipNU, "NUmagnetic_flux")

        if potential:
            _PhiNU = self.efield.PhiNU(_psiNU, theta)
            PhiNU = self.Q(_PhiNU, "NUVolts")
        else:
            PhiNU = self.Q(0, "NUVolts")

        PzetaNU = (
            2
            * self.miNU
            * gNU**2
            * (self.ENU - self.muNU * bNU - self.qiNU * PhiNU)
            / (self.qiNU * bNU) ** 2
        ) ** (1 / 2) - psipNU

        return PzetaNU.to(units)

    def findmu(
        self, psi: Quantity, theta: float, units: str, potential: bool = True
    ):

        # All Quantities and operations are in NU, and the conversion to input
        # units takes place at the end.
        psiNU = psi.to("NUMagnetic_flux")
        _psiNU = psiNU.magnitude

        # Calculate and Quantify bfield and currents
        _bNU, _iNU, _gNU = self.bfield.bigNU(_psiNU, theta)
        bNU = self.Q(_bNU, "NUTesla")
        gNU = self.Q(_gNU, "NUPlasma_current")

        # Calculate and Quantify psip
        _psipNU = self.qfactor.psipNU(_psiNU)
        psipNU = self.Q(_psipNU, "NUmagnetic_flux")

        if potential:
            _PhiNU = self.efield.PhiNU(_psiNU, theta)
            PhiNU = self.Q(_PhiNU, "NUVolts")
        else:
            PhiNU = self.Q(0, "NUVolts")

        rhoNU = (self.PzetaNU - psipNU) / gNU

        muNU = (self.ENU - self.qiNU * PhiNU) / bNU - (
            self.qiNU**2 / (2 * self.miNU)
        ) * rhoNU**2 * bNU

        return muNU.to(units)

    def __repr__(self):
        return Tokamak.__repr__(self) + PhysicalParameters.__repr__(self)

    def __str__(self):
        return Tokamak.__str__(self) + PhysicalParameters.__str__(self)
