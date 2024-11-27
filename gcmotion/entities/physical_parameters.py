r"""
=========================
PhysicalParameters Entity
=========================

This module defines the "PhysicalParameters" entity, which represents a family
of particles that are the same species, and have the same "mu" and "Pzeta"

The Physical values are stored as Quantites, as well as their NU counterparts.
"""

from __future__ import annotations

import pint
from typing import Literal
from termcolor import colored

from gcmotion.utils.logger_setup import logger
from gcmotion.configuration.physical_constants import PhysicalConstants

type SupportedSpecies = Literal["p", "e", "D", "T", "He3", "He4"]
type Quantity = pint.Quantity


class PhysicalParameters:
    r"""Returns a PhysicalParameters entity.

    Contains the constants of motion :math:`\mu` and :math:`P_\zeta`, as well
    as the particle species.

    Parameters
    ----------
    species : {'p', 'e', 'D', 'T', 'He3', 'He4'}
        The particle's species. This field is case-insensitive.
    mu : Quantity
        The magnetic moment COM in units of [current][length]^2.
    Pzeta : Quantity
        The :math:`P_\zeta` canonical momentum COM in units of magnetic
        flux.

    Notes
    -----
    A PhysicalParameters object contains the following attributes, which
    include the input arguements:

        #. mu, muNU : Quantities
            The magnetic moment in SI/NU.
        #. Pzeta, PzetaNU : Quantites
            The :math:`P_\zeta` canonical momentum COM in SI/NU.
        #. species, species_name : str
            The species abbreviated and full name.

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
    >>> Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)
    >>>
    >>> params = gcm.PhysicalParameters(
    ...     species=species,
    ...     mu=Q(1e-5, "NUMagnetic_moment"),
    ...     Pzeta=Q(-0.025, "NUMagnetic_flux"),
    ... )

    """

    def __init__(
        self, species: SupportedSpecies, mu: Quantity, Pzeta: Quantity
    ):

        logger.info("==> Initializing PhysicalParameters...")
        check(species, mu, Pzeta)

        # Construct Quantities
        self.species = species.lower()
        self.species_name = getattr(
            PhysicalConstants, self.species + "_name", None
        )
        self.mu = mu.to("Magnetic_moment")
        self.Pzeta = Pzeta.to("Magnetic_flux")

        # Corresponding NU Quantities
        self.muNU = self.mu.to("NUMagnetic_moment")
        self.PzetaNU = self.Pzeta.to("NUmagnetic_flux")

    def __repr__(self):
        return (
            f"species = {self.species}, "
            + f"mu = {self.mu:.4g~}, "
            + f"Pzeta= {self.Pzeta:.4g~}"
        )

    def __str__(self):
        return (
            colored("\nPhysical Parameters:\n", "green")
            + f"{"Particle species":>23} : "
            + f"{colored(self.species_name, "light_blue"):<16}\n"
            + f"{"mu":>23} : {f'{self.mu:.4g~}':<16}"
            + f"({self.muNU:.4g~})\n"
            + f"{"Pzeta":>23} : {f'{self.Pzeta:.4g~}':<16}"
            + f"({self.PzetaNU:.4g~})\n"
        )


def check(species, mu, Pzeta):
    r"""Checks the validity of the arguements"""

    assert species.lower() in [
        "p",
        "e",
        "D",
        "T",
        "He3",
        "He4",
    ], "species can be on of 'p', 'e', 'D', 'T', 'He3', 'He4'"
    assert isinstance(mu, pint.Quantity), "'mu' must be a Quantity!"
    assert isinstance(Pzeta, pint.Quantity), "'Pzeta' must be a Quantity"
    assert mu.dimensionality == {
        "[current]": 1,
        "[length]": 2,
    }, "'mu' must have a dimensionality of [current][length]^2!"
    assert Pzeta.dimensionality == {
        "[current]": -1,
        "[length]": 2,
        "[mass]": 1,
        "[time]": -2,
    }, (
        "'Pzeta' must have dimensionality of "
        + "{[current]^-1[length]^2[mass][time]^-2 "
        + "(Magnetic_flux)"
    )
