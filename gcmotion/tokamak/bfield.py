r"""
Magnetic Field configurations
=============================

Here is a list of the availiable bfield configurations.

==================     ===============
Large Aspect Ratio     :py:class:`LAR`
==================     ===============

Their parameters are documented below.

Example
-------

>>> import gcmotion as gcm
>>>
>>> # Quantity Constructor
>>> Rnum = 1.6
>>> anum = 0.5
>>> B0num = 1
>>> species = "p"
>>> ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)
>>>
>>> # Intermediate values
>>> B0 = Q(B0num, "Tesla")
>>> i = Q(0, "NUPlasma_current")
>>> g = Q(1, "NUPlasma_current")
>>>
>>> # A Magnetic field
>>> bfield1 = gcm.bfield.LAR(B0=B0, i=i, g=g)

The functions `bigNU` and `solverbNU` work identically in every class, so I
list their methods here as to not repeat myself:

.. autofunction:: gcmotion.bfield.MagneticField.bigNU

.. autofunction:: gcmotion.bfield.MagneticField.solverbNU

.. rubric:: LAR

.. autoclass:: LAR
    :show-inheritance:

"""

import pint
import numpy as np
from termcolor import colored

from math import cos, sin, sqrt
from abc import ABC, abstractmethod

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


class MagneticField(ABC):
    r"""Magnetic field base class"""

    def __init__(self):
        r"""Contains all the needed parameters

        The parameters should be defined as Quantities with SI units.
        """

    @abstractmethod
    def bigNU(
        self, phi: float | np.ndarray, theta: float | np.ndarray
    ) -> float | np.ndarray:
        r"""Calculates :math:`B(\psi, \theta), I(\psi, \theta), g(\psi,\
        \theta)`. Input and output must be both floats or np.ndarrays, in
        [NU].

        Used in energy contour plots.

        Parameters
        ----------
        psi : float | np.ndarray
            The :math:`\psi` value(s) in [NU].
        theta : float | np.ndarray
            The :math:`\theta` value(s).

        Returns
        -------
        float | np.ndarray
            The Calculated :math:`B(\psi, \theta), I(\psi, \theta), g(\psi,\
            \theta)` value(s) in [NU].

        """

    @abstractmethod
    def solverbNU(
        self, psi: float, theta: float
    ) -> tuple[float, float, float]:
        r"""Calculates all the values needed by the solver:
        :math:`B,I,g` (by calling ``self.bigNU()``) and the derivatives
        :math:`\dfrac{\partial B}{\partial \psi}, \dfrac{\partial B}{\partial\
        \theta}\ \dfrac{\partial I}{\partial \psi}` and :math:`\dfrac{\partial\
        g}{\partial \psi}`.

        Input and output must be floats, in [NU].

        Used inside the solver.

        .. warning::
            The derivatives are calculated with respect to :math:`\psi`, and
            **not** :math:`\psi_p`, which appear in the differential equations.
            This is accounted for inside the solver by multiplying by
            :math:`q(\psi)`, since :math:`\dfrac{\partial f}{\partial \psi_p}\
            =\dfrac{\partial f}{\partial \psi}\ \dfrac{\partial \psi}{\partial\
            \psi_p} = q\dfrac{\partial f}{\partial \psi}`

        Parameters
        ----------
        psi : float
            Value of :math:`\psi` in [NU].
        theta : float
            Value of :math:`\theta` in [NU].

        Returns
        -------
        3-tuple of floats
            Calculated values and derivatives in [NU].
        """


# =======================================================


class LAR(MagneticField):
    r"""Initializes the standard Large Aspect Ratio magnetic field.

    Parameters
    -----------
    B0 : Quantity
        The strength of the magnetic field *on the magnetic axis*
        in [T].
    i : Quantity
        The toroidal current in units of [Plasma Current].
    g : Quantity
        The poloidal current in units of [Plasma Current].

    """

    def __init__(self, B0: Quantity, i: Quantity, g: Quantity):
        r"""Parameters initialization."""

        # SI Quantities
        self.B0 = B0.to("Tesla")
        self.i = i.to("Plasma_Current")
        self.g = g.to("Plasma_Current")

        # [NU] Conversions
        self.B0NU = B0.to("NUTesla")
        self.iNU = i.to("NUPlasma_Current")
        self.gNU = g.to("NUPlasma_Current")

        # Unitless quantities, makes it a bit faster if defined here
        for key, value in self.__dict__.copy().items():
            self.__setattr__("_" + key, value.magnitude)

        # Flags
        self.has_i = bool(i.m)

    def bigNU(self, psi: float | np.ndarray, theta: float | np.ndarray):

        if isinstance(psi, (int, float)):
            b = 1 - sqrt(2 * psi) * cos(theta)
            g = self._gNU
            i = self._iNU
        elif isinstance(psi, np.ndarray):
            b = 1 - np.sqrt(2 * psi) * np.cos(theta)
            g = np.tile(self._gNU, psi.shape)
            i = np.tile(self._iNU, psi.shape)

        return (b, i, g)

    def solverbNU(self, psi: float, theta: float):

        # Field and currents
        b, i, g = self.bigNU(psi, theta)

        # Field derivatives
        root = sqrt(2 * psi)
        cos_theta = cos(theta)
        b_der_psi = cos_theta / root
        b_der_theta = root * sin(theta)

        # Current derivatives
        i_der = 0
        g_der = 0

        # Pack them up
        currents = (i, g)
        b_der = (b_der_psi, b_der_theta)
        currents_der = (i_der, g_der)

        return b, b_der, currents, currents_der

    def __repr__(self):
        return (
            colored("LAR", "light_blue")
            + f": B0={self.B0:.4g~}, I={self.i:.4g~}, g={self.g:.4g~}."
        )
