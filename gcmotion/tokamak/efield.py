"""
About Electric field objects
============================

An Electric Field object is a class instance containing all the information 
about the electric field of the system. It implements all the methods needed 
buy the solver and other calculations, and is called automatically wherever required.

To add a new Electric Field, simply copy-paste an already existing class
(idealy the Nofield one) and fill the ``__init__()``, ``Phi_der()``,
``Er_of_psi()``, and ``Phi_of_psi()`` methods to fit your Electric Field.  
In case your Electric field has extra parameters you want to pass as 
arguments, you must also create an ``__init__()`` method and declare them. 
To avoid errors, your class should inherit the ``ElectricField`` class.

.. caution::
    **All values, both input and output are in SI units.**
    
    Specifically [V/m] and [V].

    The conversion is done internally.

The general structure is this::

    class MyElectricField(ElectricField):

        def __init__(self, *<parameters>):
            self.id = "foo" # Simple id used only for logging.
            self.params = {} # Tweakable parameters, used only for logging.
            <set parameters>
    
        def Phi_der(self, psi): 
            return [Phi_der_psip, Phi_der_theta]
    
        def Er_of_psi(self, psi):
            r = np.sqrt(2 * psi)
            return E
    
        def Phi_of_psi(self, psi):
            return Phi

.. note::
    The above methods return the same type as the input (either Python floats 
    or np.ndarrays). When used inside the solver, they should return a 
    Python float, and not a np.float. Solvers need to be fast so they work 
    with built-in floats, while plotting functions work with np.ndarrays.This 
    is mainly for optimization reasons and should probably not cause problems.

.. rubric:: The 'ElectricField' Abstract Base Class

The base class that every other class inherits from is ``ElectricField``. 
This class does nothing, it is only a template.

.. autoclass:: ElectricField
    :members: __init__, Phi_der, Er_of_psi, Phi_of_psi

:show-inheritance:
"""

import numpy as np
from scipy.special import erf
from .qfactor import QFactor
from math import sqrt, exp
from abc import ABC, abstractmethod


class ElectricField(ABC):
    r"""Electric field base class."""

    def __init__(self):
        self.id = "Base Class"
        self.params = {}

    @abstractmethod
    def Phi_der(self, psi: float) -> tuple[float, float]:
        r"""Calculates the derivatives of Φ(ψ) with respect to
        :math:`\psi_p, \theta` in [V].

        Intended for use only inside the ODE solver. Returns the potential
        in [V], so the normalisation is done inside the solver.

        Parameters
        ----------

        psi : float
            The magnetic flux surface.

        Returns
        -------

        tuple : 2-tuple of floats
            2-tuple containing the calculated derivatives.
        """
        pass

    @abstractmethod
    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r"""Calculates radial Electric field component in [V/m] from ψ.

        Used for plotting the Electric field

        Parameters
        ----------

        psi : np.ndarray
            The ψ values.

        Returns
        -------

        E : np.ndarray
            1D numpy array with calculated E values.

        """
        pass

    @abstractmethod
    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r"""Calculates Electric Potential in [V] from ψ.

        Used for plotting the Electric Potential, the particles initial Φ,
        and the Φ values for the contour plot.

        Parameters
        ----------

        psi : np.ndarray
            The ψ values.

        Returns
        -------

        Phi : np.ndarray
            1D numpy array with calculated values.

        """
        pass


# =======================================================


class Nofield(ElectricField):
    r"""Initializes an electric field of 0

    Exists to avoid compatibility issues.

    Takes no parameters.
    """

    def __init__(self):
        self.id = "NoField"
        self.params = {}
        return

    def Phi_der(self, psi: float) -> tuple[float, float]:
        return (0, 0)

    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        return 0 * psi

    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        return 0 * psi


class Parabolic(ElectricField):
    r"""Initializes an electric field of the form: :math:`E(r) = ar^2 + b`"""

    def __init__(self, R: float, a: float, q: QFactor, alpha: float, beta: float):
        r"""Initializes the field's parameters.

        Parameters
        ----------

        R : float
            The tokamak's major radius in [m].
        a : float
            The tokamak's minor radius in [m].
        q : :class:`QFactor` object
            q factor profile.
        alpha : float
            The :math:`r^2` coefficient.
        beta : float
            The constant coefficient.

        """
        self.id = "Parabolic"
        self.params = {"alpha": alpha, "beta": beta}

        self.a = alpha
        self.b = beta
        self.q = q
        self.r_wall = a / R
        self.psi_wall = (self.r_wall) ** 2 / 2  # normalized to R
        self.psip_wall = q.psip_of_psi(self.psi_wall)

    def Phi_der(self, psi: float) -> tuple[float, float]:
        r = np.sqrt(2 * psi)
        Phi_der_psip = -self.q.q_of_psi(psi) * (self.a * r - self.b / r)
        Phi_der_theta = 0
        return [Phi_der_psip, Phi_der_theta]

    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r = np.sqrt(2 * psi)
        E = self.a * r**2 + self.b
        return E

    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r = np.sqrt(2 * psi)
        Phi = -(self.a / 3) * r**3 - self.b * r
        return Phi


class Radial(ElectricField):
    r"""Initializes an electric field of the form:
    :math:`E(r) = -E_a\exp\bigg[-\dfrac{(r-r_a)^2}{r_w^2})\bigg]`
    """

    def __init__(
        self,
        R: float,
        a: float,
        q: QFactor,
        Ea: float,
        minimum: float,
        waist_width: float,
    ):
        r"""Initializes the field's parameters.

        Parameters
        ----------

        R : int | float
            The tokamak's major radius in [m].
        a : int | float
            The tokamak's minor radius in [m].
        q : :class:`QFactor` object
            q factor profile.
        Ea : float
            The Electric field magnitude in [V/m].
        minimum : float
            The Electric field's minimum point with respect to
            :math:`\psi_{wall}`.
        waist_width : float
            The Electric field's waist width, defined as:
            :math:`r_w = \dfrac{a}{\text{waste width}}`.

        """
        self.id = "Radial"
        self.params = {"Ea": Ea, "minimum": minimum, "waist_width": waist_width}

        self.q = q
        self.r_wall = a / R
        self.psi_wall = (self.r_wall) ** 2 / 2  # normalized to R
        self.psip_wall = q.psip_of_psi(self.psi_wall)
        self.minimum = minimum
        self.waist_width = waist_width

        self.Ea = Ea  # V/m
        self.ra = self.minimum * self.r_wall  # Defines the minimum point
        self.Efield_min = self.ra**2 / 2
        self.rw = self.r_wall / self.waist_width  # waist, not wall
        self.psia = self.ra**2 / 2
        self.psiw = self.rw**2 / 2  # waist, not wall

        # Square roots, makes it a bit faster
        self.sr_psia = sqrt(self.psia)
        self.sr_psiw = sqrt(self.psiw)

    def Phi_der(self, psi: float) -> tuple[float, float]:
        Phi_der_psip = (
            self.q.q_of_psi(psi)
            * self.Ea
            / (sqrt(2 * psi))
            * exp(-((sqrt(psi) - self.sr_psia) ** 2) / self.psiw)
        )
        Phi_der_theta = 0

        return [Phi_der_psip, Phi_der_theta]

    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r = np.sqrt(2 * psi)
        E = -self.Ea * np.exp(-((r - self.ra) ** 2) / self.rw**2)
        return E

    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        Phi = (
            self.Ea
            * np.sqrt(np.pi * self.psiw / 2)
            * (erf((np.sqrt(psi) - self.sr_psia) / self.sr_psiw) + erf(self.sr_psia / self.sr_psiw))
        )
        return Phi
