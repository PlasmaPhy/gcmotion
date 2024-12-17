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
from termcolor import colored
from typing import Literal

from gcmotion.utils.logger_setup import logger

from gcmotion.entities.tokamak import Tokamak
from gcmotion.configuration.physical_constants import PhysicalConstants

type Quantity = pint.Quantity
type SupportedSpecies = Literal["p", "e", "D", "T", "He3", "He4"]


class Profile:
    r"""A Profile entity describes an unperturbed equilibrium with a given
    Tokamak configuration and with at least 2/3 fixed Constants of Motion.

    By fixing 2 COMs and letting the other one be variable, we can perform many
    useful analyses.

    Note that even if specified, the 3rd COM will be ignored in an analysis
    that by its nature lets it vary. For example, contouring over Energy
    levels, will simply ignore the E parameter.


    Parameters
    ----------
    tokamak : Tokamak
        The Tokamak entity.
    species : str
        The particles' species.
    mu : Quantity
        The Magnetic Moment COM :math:`\mu`.
    Pzeta : Quantity
        The Canonical Momentum COM :math:`P_\zeta`.
    E : Quantity
        The E COM :math:`E`.

    Example
    -------
    How to create a `Profile` object.

    >>> import gcmotion as gcm
    >>>
    >>> # Quantity Constructor
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
    >>> mu = Q(1e-5, "NUMagnetic_moment")
    >>> Pzeta = Q(-0.015, "NUMagnetic_flux")
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
    >>> # Create a Profile
    >>> profile = gcm.Profile(
    ...     tokamak=tokamak,
    ...     species=species,
    ...     mu=mu,
    ...     Pzeta=Pzeta,
    ... )

    """

    def __init__(
        self,
        tokamak: Tokamak,
        species: SupportedSpecies,
        mu: Quantity = None,
        Pzeta: Quantity = None,
        E: Quantity = None,
    ):

        logger.info("==> Initializing Profile...")

        # Check if at least 2 are given
        if [mu, Pzeta, E].count(None) > 1:
            msg = "At least 2/3 Constants of motion must be specified"
            raise ValueError(msg)

        # Grab attributes from tokamak object
        self.tokamak = tokamak
        self.__dict__.update(self.tokamak.__dict__)

        # Define species
        self.Q = type(tokamak.R)
        self.species = species.lower()
        self.species_name = getattr(
            PhysicalConstants, self.species + "_name", None
        )

        # Grab particle's mass and charge
        M = getattr(PhysicalConstants, self.species + "_M")
        Z = getattr(PhysicalConstants, self.species + "_Z")

        self.miNU = self.Q(M, "Proton_mass")
        self.qiNU = self.Q(Z, "Proton_charge")

        self.mi = self.miNU.to("kilogram")
        self.qi = self.qiNU.to("Coulomb")

        # Constants of motion
        if mu is None:
            self.mu = self.muNU = None
        else:
            self.mu = mu.to("Magnetic_moment")
            self.muNU = self.mu.to("NUMagnetic_moment")

        if Pzeta is None:
            self.Pzeta = self.PzetaNU = None
        else:
            self.Pzeta = Pzeta.to("Magnetic_flux")
            self.PzetaNU = self.Pzeta.to("NUmagnetic_flux")

        if E is None:
            self.E = self.ENU = None
        else:
            self.E = E.to("keV")
            self.ENU = self.E.to("NUJoule")

        self.Q = type(self.tokamak.R)

        logger.info("\tProfile setup complete.")

    def findPtheta(self, psi: Quantity):
        r"""Calculates Ptheta from psi. ``Pzeta`` must be defined.

        Output units are the same as input units.

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
        pair. Both `Pzeta` and `mu` must be defined.

        Parameters
        ----------
        psi : Quantity
            The particle's psi Quantity.
        theta : float
            The particle's :math:`\theta` angle.
        units : str
            The returned Energy units.
        potential : bool, optional
            Whether or not to add the electric potential term in the energy.
            Defaults to True.

        Returns
        -------
        Quantity
            The calculated Energy Quantity in the specified units.
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
        r"""Calculates the :math:`P_\zeta` COM. Both "mu" and "E" must be
        defined.

        Parameters
        ----------
        psi : Quantity
            The particle's psi Quantity.
        theta : float
            The particle's :math:`\theta` angle.
        units : str
            The returned Pzeta units.
        potential : bool, optional
            Whether or not to add the electric potential term in the
            calculation. Defaults to True.

        Returns
        -------
        Quantity
            The calculated Pzeta Quantity in the specified units.

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
        r"""Calculates the :math:`\mu` COM. Both "Pzeta" and "E" must be
        defined.

        Parameters
        ----------
        psi : Quantity
            The particle's psi Quantity.
        theta : float
            The particle's :math:`\theta` angle.
        units : str
            The returned mu units.
        potential : bool, optional
            Whether or not to add the electric potential term in the
            calculation. Defaults to True.

        Returns
        -------
        Quantity
            The calculated mu Quantity in the specified units.
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
        string = Tokamak.__repr__(self)
        string += (
            "PhysicalParameters: "
            + f"species = {self.species}, "
            + f"mu = {self.mu}, "
            + f"Pzeta = {self.Pzeta}, "
            + f"E = {self.E}\n"
        )
        return string

    def __str__(self):
        mu = "None" if self.mu is None else f"{self.mu:.4g~}"
        muNU = "None" if self.mu is None else f"{self.muNU:.4g~}"
        Pzeta = "None" if self.Pzeta is None else f"{self.Pzeta:.4g~}"
        PzetaNU = "None" if self.Pzeta is None else f"{self.PzetaNU:.4g~}"
        E = "None" if self.E is None else f"{self.E:.4g~}"
        ENU = "None" if self.ENU is None else f"{self.ENU:.4g~}"

        string = Tokamak.__str__(self)
        string += (
            colored("\nConstans of Motion\n", "green")
            + f"{"Particle species":>23} : "
            + f"{colored(self.species_name, "light_blue"):<16}\n"
            + f"{"mu":>23} : {mu:<16}({muNU})\n"
            + f"{"Pzeta":>23} : {f'{Pzeta}':<16}({PzetaNU})\n"
            + f"{"E":>23} : {f'{E}':<16}({ENU})\n"
        )
        return string
