r"""
q-factor configurations
=======================

Here is a list of the availiable q-factor configurations:

=======================      ==========================
Unity q-factor               :py:class:`Unity`
Parabolic q-factor           :py:class:`Parabolic`
Hypergeometric q-factor      :py:class:`Hypergeometric`
=======================      ==========================

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
>>> a = Q(anum, "meters")
>>> B0 = Q(B0num, "Tesla")
>>>
>>> # Some qfactors
>>> qfactor1 = gcm.qfactor.Unity()
>>> qfactor3 = gcm.qfactor.Parabolic(a, B0, q0=1.1, q_wall=3.8)
>>> qfactor3 = gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2)

.. note::

    The values of `a` and `B0` are necessary whenever we define the qfactor
    with reference to its value at the wall.

The functions `solverqNU` and `psipNU` work identically in every class, so I
list their methods here as to not repeat myself:

.. autofunction:: gcmotion.qfactor.QFactor.solverqNU

.. autofunction:: gcmotion.qfactor.QFactor.psipNU


.. rubric:: Unity

.. autoclass:: Unity
    :member-order: bysource
    :inherited-members: solverqNU, psipNU
    :undoc-members: solverqNU, psipNU
    :show-inheritance:

.. rubric:: Parabolic

.. autoclass:: Parabolic
    :show-inheritance:

.. rubric:: Hypergeometric

.. autoclass:: Hypergeometric
    :show-inheritance:

"""

import pint
import numpy as np

from termcolor import colored
from math import sqrt, atan
from scipy.special import hyp2f1
from abc import ABC, abstractmethod

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


class QFactor(ABC):
    r"""q Factor base class"""

    @abstractmethod
    def __init__(self):
        r"""Contains all the needed parameters."""

    @abstractmethod
    def solverqNU(self, psi: float) -> float:
        r"""Calculates :math:`q(\psi)`.
        Input and output must both be floats, in [NU].

        Used inside the solver.

        Parameters
        ----------
        psi : float
            Value of :math:`\psi` in [NU].

        Returns
        -------
        float
            Calculated :math:`q(\psi)`.

        """

    @abstractmethod
    def psipNU(self, psi: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates :math:`\psi_p(\psi)`.
        Input and output must can be both floats or np.ndarrays, in [NU].

        Used in calculating :math:`\psi_p`'s time evolution from :math:`\psi`,
        :math:`\psi_{p,wall}`, etc.

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
    r"""Initializes an object q with :math:`q(\psi) = 1`
    and :math:`\psi_p=\psi`."""

    def __init__(self):
        pass

    def solverqNU(self, psi):
        r"""Always returns 1."""
        return 1

    def psipNU(self, psi):
        """Always returns `psi`."""
        return psi

    def __repr__(self):
        return colored("Unity", "light_blue")


class Parabolic(QFactor):
    r"""Initializes an object q with

    .. math::

        q(\psi) = q_0 + (q_{wall}-q_0)\bigg(\dfrac{\psi}{\psi_{wall}}\
        \bigg)^2

    :math:`\psi_p(\psi)` is calculated from:

    .. math::

        \psi_p(\psi) = \dfrac{\psi_{wall}}{\sqrt{q_0 (q_{wall}-q_0)}}\
        \arctan\bigg( \dfrac{\psi\sqrt{q_{wall}-q_0}}{\psi_{wall}\sqrt{q_0}}\
        \bigg)

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

    def __init__(
        self,
        a: Quantity,
        B0: Quantity,
        q0: float,
        q_wall: float,
    ):
        r"""Parameters initialization."""

        # SI Quantities
        self.B0 = B0.to("Tesla")
        self.a = a.to("meters")
        self.psi_wall = (B0 * a**2 / 2).to("Magnetic_flux")
        self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")

        # Unitless quantities, makes it a bit faster if defined here
        self._psi_wall = self.psi_wall.magnitude
        self._psi_wallNU = self.psi_wallNU.magnitude
        self.sra = sqrt(q0)
        self.srb = sqrt(q_wall - q0)

        # Purely Numerical Parameters
        self.q0 = q0
        self.q_wall = q_wall

    def solverqNU(self, psi):
        return (
            self.q0 + (self.q_wall - self.q0) * (psi / self._psi_wallNU) ** 2
        )

    def psipNU(self, psi):
        if isinstance(psi, float):
            return (self._psi_wallNU / (self.sra * self.srb)) * atan(
                self.srb * psi / (self.sra * self._psi_wallNU)
            )
        else:
            return (self._psi_wallNU / (self.sra * self.srb)) * np.atan(
                self.srb * psi / (self.sra * self._psi_wallNU)
            )

    def __repr__(self):
        return (
            colored("Parabolic", "light_blue")
            + f": q0={self.q0:.4g}, q_wall={self.q_wall:.4g}."
        )


class Hypergeometric(QFactor):
    r"""Initializes an object q with:

    .. math::

        q(\psi) = q_0\bigg\{ 1 + \bigg[ \bigg(\dfrac{q_{wall}}{q_0}\
        \bigg)^n -1 \bigg] \
        \bigg( \dfrac{\psi}{\psi_{wall}} \bigg)^n \bigg\}^{1/n}

    :math:`\psi_p(\psi)` is calculated from:

    .. math::

        \psi_p(\psi) = \dfrac{\psi}{q_0} \phantom{1}_2 F_1\
        \bigg[ \dfrac{1}{n}, \dfrac{1}{n}, 1+\dfrac{1}{n},
        \bigg(1 - \bigg( \dfrac{q_{wall}}{q_0} \bigg)^n\bigg)
        \bigg( \dfrac{\psi}{\psi_{wall}} \bigg)^n \bigg]

    where :math:`\phantom{1}_2 F_1` the hypergeometric function.

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

        # [NU] Conversions
        self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")

        # Unitless quantities, makes it a bit faster if defined here
        for key, value in self.__dict__.copy().items():
            self.__setattr__("_" + key, value.magnitude)

        # Purely Numerical Quantities
        self.q0 = q0
        self.q_wall = q_wall
        self.n = n

    def solverqNU(self, psi):
        return self.q0 * (
            1
            + ((self.q_wall / self.q0) ** self.n - 1)
            * (psi / self._psi_wallNU) ** self.n
        ) ** (1 / self.n)

    def psipNU(self, psi):

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
            colored("Hypergeometric", "light_blue")
            + f": q0={self.q0:.4g}, q_wall={self.q_wall:.4g}, n={self.n:.4g}."
        )
