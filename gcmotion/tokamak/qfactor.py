r"""
About q-Factor objects
----------------------

A q-Factor object is a class instance containing all the information about the 
q-factor profile of the system. It implements all the methods needed buy the
solver and other calculations, and is called automatically wherever required.

To add a new q-factor, simply copy-paste an already existing class
(idealy the Unity one) and fill the ``q_of_psi()`` and ``psip_of_psi()`` 
methods to fit your q-factor. In case your q factor has extra parameters
you want to pass as arguments, you must also create an ``__init__()``
method and declare them. To avoid errors, your class should inherit the
``QFactor`` class.

The general structure is this::

    class MyQFactor(QFactor):

        def __init__(self, *parameters):
            "Parameter setup."

        def q_of_psi(self, psi):
            "Returns the value q(ψ)."
            return q

        def psip_of_psi(self, psi):
            "Returns the value ψ_p(ψ)."
            return pisp

.. note::
    The Qfactor's attributes must be Quantities with units in SI so they can be
    referenced from anywhere in the code, but the methods input **must** be purely
    numeric. As a result, the output is also purely numeric and in the same units 
    as the input.

.. important::
    q-factor is a dimensionless quantity, and values such as 
    :math:`q_0, q_{wall}, \psi, \psi_{wall}, \ldots` *seem* to always appear in 
    ratios. Therefore we don't need to convert anything to [NU] when calculating
    q inside the solver, since the result will be the same!

.. rubric:: The 'QFactor' Abstract Base Class

The base class that every other class inherits from is ``QFactor``. This class 
does nothing, it is only a template.

.. autoclass:: QFactor
    :member-order: bysource
    :members: __init__, q_of_psi, psip_of_psi

"""

import pint
import numpy as np
from math import sqrt, atan
from scipy.special import hyp2f1
from abc import ABC, abstractmethod

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


class QFactor(ABC):
    r"""q Factor base class"""

    @abstractmethod
    def __init__(self):
        r"""Contains all the needed parameters

        The parameters should be defined as Quantities with SI units.
        """

    @abstractmethod
    def q_of_psi(self, psi: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates q(ψ). Output size is the same as input size.

        When used inside the solver, it returns :math:`q(\psi)` as a float
        When used for the tokamak profile plotting, it returns an np.ndarray.

        Parameters
        ----------
        psi : float | np.ndarray
            Value(s) of ψ.

        Returns
        -------
        float | np.ndarray
            Calculated q(ψ)

        """

    @abstractmethod
    def psip_of_psi(self, psi: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates :math:`\psi_p(\psi)`. Output size is the same as input size.

        Used in calculating :math:`\psi_{p,wall}`, :math:`\psi_p`'s
        time evolution from :math:`\psi`, etc.

        Parameters
        ----------
        psi : float | np.ndarray
            Value(s) of ψ.

        Returns
        -------
        float | np.ndarray
            Calculated :math:`\psi_p(\psi)`.

        """


# ====================================================


class Unity(QFactor):
    r"""Initializes an object q with :math:`q(\psi) = 1`"""

    def __init__(self):
        pass

    def q_of_psi(self, psi):
        return 1

    def q_of_psiNU(self, psi):
        return 1

    def psip_of_psi(self, psi):
        return psi

    def __repr__(self):
        return "Unity"


class Parabolic(QFactor):
    r"""Initializes an object q with
    :math:`q(\psi) = q_0 + (q_{wall}-q_0)\bigg(\dfrac{\psi}{\psi_{wall}}\bigg)^2`
    """

    def __init__(
        self,
        a: Quantity,
        B0: Quantity,
        q0: float,
        q_wall: float,
    ):
        r"""Parameters initialization.

        Parameters
        ----------
        a : Quantity
            float The tokamak's minor radius in [m].
        B0 : Quantity
            The Magnetic field's strength in [T]
        q0 : float
            q-value at the magnetic axis.
        q_wall : float
            q_value at the wall.
        n : int
            Order of equillibrium (1: peaked, 2: round, 4: flat).
        """

        # SI Quantities
        self.B0 = B0.to("Tesla")
        self.a = a.to("meters")
        self.psi_wall = (B0 * a**2 / 2).to("Magnetic_flux")

        # [NU] Conversions for q_of_psiNU():
        self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")

        # Unitless quantities, makes it a bit faster if defined here
        self._psi_wall = self.psi_wall.magnitude
        self._psi_wallNU = self.psi_wallNU.magnitude
        self.sra = sqrt(q0)
        self.srb = sqrt(q_wall - q0)

        # Purely Numerical Parameters
        self.q0 = q0
        self.q_wall = q_wall

    def q_of_psi(self, psi):
        return (
            self.q0 + (self.q_wall - self.q0) * (psi / self._psi_wallNU) ** 2
        )

    def psip_of_psi(self, psi):
        if isinstance(psi, float):
            return (self._psi_wallNU / (self.sra * self.srb)) * atan(
                self.srb * psi / (self.sra * self._psi_wallNU)
            )
        else:
            return (self._psi_wallNU / (self.sra * self.srb)) * np.atan(
                self.srb * psi / (self.sra * self._psi_wallNU)
            )

    def __repr__(self):
        return "Parabolic: " + f"q0={self.q0:.4g}, q_wall={self.q_wall:.4g}."


class Hypergeometric(QFactor):
    r"""Initializes an object q with
    :math:`q(\psi) = q_0\bigg\{ 1 + \bigg[ \bigg(\dfrac{q_{wall}}{q_0}\bigg)^n -1 \bigg] \
    \bigg( \dfrac{\psi}{\psi_{wall}} \bigg)^n \bigg\}^{1/n}`.
    """

    def __init__(
        self,
        a: Quantity,
        B0: Quantity,
        q0: float,
        q_wall: float,
        n: int,
    ):
        r"""Parameters initialization.

        Parameters
        ----------
        a : Quantity
            float The tokamak's minor radius in [m].
        B0 : Quantity
            The Magnetic field's strength in [T].
        q0 : float
            q-value at the magnetic axis.
        q_wall : float
            q_value at the wall.
        n : int
            Order of equillibrium (1: peaked, 2: round, 4: flat).
        """

        # SI Quantities
        self.B0 = B0.to("Tesla")
        self.a = a.to("meters")
        self.psi_wall = (B0 * a**2 / 2).to("Magnetic_flux")

        # [NU] Conversions for q_of_psiNU():
        self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")

        # Unitless quantities, makes it a bit faster if defined here
        self._psi_wall = self.psi_wall.magnitude
        self._psi_wallNU = self.psi_wallNU.magnitude

        # Purely Numerical Quantities
        self.q0 = q0
        self.q_wall = q_wall
        self.n = n

    def q_of_psi(self, psi):
        return self.q0 * (
            1
            + ((self.q_wall / self.q0) ** self.n - 1)
            * (psi / self._psi_wallNU) ** self.n
        ) ** (1 / self.n)

    def psip_of_psi(self, psi):

        a = b = 1 / self.n
        c = 1 + 1 / self.n
        z = (1 - (self.q_wall / self.q0) ** self.n) * (
            psi / self._psi_wallNU
        ) ** self.n
        if isinstance(psi, (int, float)):
            return psi / self.q0 * float(hyp2f1(a, b, c, z))
        else:
            return psi / self.q0 * hyp2f1(a, b, c, z)

    def __repr__(self):
        return (
            "Hypergeometric: "
            + f"q0={self.q0:.4g}, q_wall={self.q_wall:.4g}, n={self.n:.4g}."
        )
