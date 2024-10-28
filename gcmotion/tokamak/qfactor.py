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

        def __init__(self, *<parameters>):
            self.id = "foo" # Simple id used only for logging.
            self.params = {} # Tweakable parameters, used only for logging.
            <set parameters>

        def q_of_psi(self, psi):
            "Returns the value q(ψ)."
            return q

        def psip_of_psi(self, psi):
            "Returns the value ψ_p(ψ)."
            return pisp

.. note::
    Keep in mind that when these two methods return the same type as the input 
    (either Python floats or np.ndarrays). When used inside the solver, they 
    should return a Python float, and not a np.float. Solvers need to be fast 
    so they work with built-in floats, while plotting functions work with 
    np.ndarrays.This is mainly for optimization reasons and should probably 
    not cause problems.

.. rubric:: The 'QFactor' Abstract Base Class

The base class that every other class inherits from is ``QFactor``. This class 
does nothing, it is only a template.

.. autoclass:: QFactor
    :members: __init__, q_of_psi, psip_of_psi

"""

import numpy as np
from math import atan
from scipy.special import hyp2f1
from abc import ABC, abstractmethod


class QFactor(ABC):
    r"""q Factor base class"""

    @abstractmethod
    def __init__(self):
        r"""Contains all the needed parameters"""
        self.id = "Base Class"
        self.params = {}

    @abstractmethod
    def q_of_psi(self, psi: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates q(ψ). Return type should be same as input.

        Used inside dSdt, Φ derivatives (returns a float) and plotting of
        q factor (returns an np.ndarray).

        Parameters
        ----------

        psi : float | np.ndarray
            Value(s) of ψ.

        Returns
        -------

        qfactor : float | np.ndarray
            Calculated q(ψ)

        """
        pass

    @abstractmethod
    def psip_of_psi(self, psi: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates :math:`\psi_p(\psi)`.

        Used in calculating :math:`\psi_{p,wall}` in many methods (returns a float), in calculating
        :math:`\psi_p`'s time evolution (returns an np.ndarray), in Energy contour calculation
        (returns an np.ndarray) and in q-factor plotting (returns an np.ndarray)

        Parameters
        ----------

        psi : float | np.ndarray
            Value(s) of ψ.

        Returns
        -------

        psip : float | np.ndarray
            Calculated :math:`\psi_p(\psi)`.

        """
        pass


# ====================================================


class Unity(QFactor):
    r"""Initializes an object q with :math:`q(\psi) = 1`"""

    def __init__(self):
        self.id = "Unity"
        self.params = {}
        return

    def q_of_psi(self, psi):
        return 1

    def psip_of_psi(self, psi):
        return psi


class Parabolic(QFactor):
    r"""Initializes an object q with :math:`q(\psi) = 1 + \psi^2`"""

    def __init__(self):
        self.id = "Parabolic"
        self.params = {}
        return

    def q_of_psi(self, psi):
        return 1 + psi**2

    def psip_of_psi(self, psi):
        return atan(psi)


class Hypergeometric(QFactor):
    r"""Initializes an object q with
    :math:`q(\psi) = q_0\bigg[ 1 + \bigg( \dfrac{\psi}{\psi_k(q_{wall})} \bigg)^n \bigg]^{1/n}`.
    """

    def __init__(
        self,
        R: int | float,
        a: int | float,
        q0: int | float,
        q_wall: int | float,
        n: int,
    ):
        r"""Parameters initialization.

        Parameters
        ----------

        R : int, float
            The tokamak's major radius.
        a : int, float
            float The tokamak's minor radius.
        q0 : int, float
            q-value at the magnetic axis.
        q_wall : int, float
            q_value at the wall
        n : int
            Order of equillibrium (1: peaked, 2: round, 4: flat).
        """
        self.id = "Hypergeometric"
        self.params = {"q0": q0, "q_wall": q_wall, "n": n}

        self.r_wall = a / R  # normalized to R
        self.psi_wall = (self.r_wall) ** 2 / 2
        self.q0 = q0
        self.q_wall = q_wall
        self.psi_knee = 0.28 * self.psi_wall
        self.n = n

    def q_of_psi(self, psi):
        # return self.q0 * (1 + (psi / (self.psi_knee)) ** self.n) ** (
        #     1 / self.n
        # )
        return self.q0 * (
            1
            + ((self.q_wall / self.q0) ** self.n - 1)
            * (psi / self.psi_wall) ** self.n
        ) ** (1 / self.n)

    def psip_of_psi(self, psi):
        a = b = 1 / self.n
        c = 1 + 1 / self.n
        z = (1 - (self.q_wall / self.q0) ** self.n) * (
            psi / self.psi_wall
        ) ** self.n
        if type(psi) is float:
            return psi / self.q0 * float(hyp2f1(a, b, c, z))
        else:
            return psi / self.q0 * hyp2f1(a, b, c, z)
