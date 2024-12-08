r"""
=========================
PhysicalParameters Entity
=========================

This module defines the "PhysicalParameters" entity, which represents a family
of particles that are the same species, and have the same "mu", "Pzeta" and
"Energy".

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
    r"""Creates a set specifying a configuration of the 3 Constants of Motion,
    :math:`\mu, P_\zeta` and :math:`E`, as well as a particle species.

    Some plots or analysis require 2 fixed COMs and one varying. In that case,
    the 3rd COM is simply ignored. For example, plotting the Energy contour of
    the profile trough
    :py:func:`~gcmotion.plot.profile_Energy_contour`, the object's *E*
    attribute is ignored.

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
        self,
        species: SupportedSpecies,
        mu: Quantity = None,
        Pzeta: Quantity = None,
        E: Quantity = None,
    ):

        logger.info("==> Initializing PhysicalParameters...")

        # Check if at only 2 are given
        if [mu, Pzeta, E].count(None) > 1:
            msg = "At least 2/3 Constants of motion must be specified"
            raise ValueError(msg)

        # Define species
        self.species = species.lower()
        self.species_name = getattr(
            PhysicalConstants, self.species + "_name", None
        )

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

    def __repr__(self):
        return (
            "PhysicalParameters: "
            + f"species = {colored(self.species, "light_blue")}, "
            + f"mu = {self.mu}, "
            + f"Pzeta = {self.Pzeta}, "
            + f"E = {self.E}\n"
        )

    def __str__(self):
        mu = "None" if self.mu is None else f"{self.mu:.4g~}"
        muNU = "None" if self.mu is None else f"{self.muNU:.4g~}"
        Pzeta = "None" if self.Pzeta is None else f"{self.Pzeta:.4g~}"
        PzetaNU = "None" if self.Pzeta is None else f"{self.PzetaNU:.4g~}"
        E = "None" if self.E is None else f"{self.E:.4g~}"
        ENU = "None" if self.ENU is None else f"{self.ENU:.4g~}"

        return (
            colored("\nPhysical Parameters:\n", "green")
            + f"{"Particle species":>23} : "
            + f"{colored(self.species_name, "light_blue"):<16}\n"
            + f"{"mu":>23} : {mu:<16}({muNU})\n"
            + f"{"Pzeta":>23} : {f'{Pzeta}':<16}({PzetaNU})\n"
            + f"{"E":>23} : {f'{E}':<16}({ENU})\n"
        )
