"""
About Electric field objects
----------------------------

An Electric Field object is a class instance containing all the information 
about the electric field of the system. It implements all the methods needed 
buy the solver and other calculations, and is called automatically wherever required.

To add a new Electric Field, simply copy-paste an already existing class
(idealy the Nofield one) and fill the ``__init__()``, ``Phi_der()``,
``Er_of_psi()``, and ``Phi_of_psi()`` methods to fit your Electric Field.  
In case your Electric field has extra parameters you want to pass as 
arguments, you must also create an ``__init__()`` method and declare them. 
To avoid errors, your class should inherit the ``ElectricField`` class.

The general structure is this::

    class MyElectricField(ElectricField):

        def __init__(self, *parameters):
            "Parameter setup."
    
        def Phi_der(self, psi): 
            return [Phi_der_psip, Phi_der_theta]
    
        def Er_of_psi(self, psi):
            return E
    
        def Phi_of_psi(self, psi):
            return Phi

.. note::
    The Electric Fields's attributes must be Quantities with units in SI 
    so they can be referenced from anywhere in the code, but the methods
    input **must** be purely numeric. As a result, the output is also 
    purely numeric and in the same units as the input.

.. caution::
    The :py:meth:`~gcmotion.tokamak.efield.ElectricField.Phi_der` method
    uses the attributes **converted to [NU]**, since it is **only** used
    inside the solver.

.. admonition:: For developers

    For each attribute that is defined as a Quantity with SI units, another,
    "hidden" attribite is also defined as its magnitude. This hidden attribute
    is then used for all the purely numerical calculations. For example:

    .. code-block:: python

        def __init__(...):
            self.rpeak = self.a * self.peak
            ...
            self._rpeak = self.rpeak.magnitude

    Here, :code:`self.rpeak` is a Quantity with units of "meters", however only
    :code:`self._rpeak` is used inside the methods. Also, by defining it in 
    :code:`__init__()` we avoid having to retrieve its magnitude every time
    a method that needs it is called.

.. rubric:: The 'ElectricField' Abstract Base Class

The base class that every other class inherits from is ``ElectricField``. 
This class does nothing, it is only a template.

.. autoclass:: ElectricField
    :member-order: bysource
    :members: __init__, Phi_der, Er_of_psi, Phi_of_psi

"""

import pint
import numpy as np
from scipy.special import erf
from math import sqrt, exp
from abc import ABC, abstractmethod

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


class ElectricField(ABC):
    r"""Electric field base class."""

    def __init__(self):
        r"""Contains all the needed parameters

        The parameters should be defined as Quantities with SI units.
        """

    @abstractmethod
    def e_values(self, psi: float) -> tuple[float, float]:
        r"""Derivatives of :math:`\Phi(\psi)` with respect to
        :math:`r, \theta`.

        Intended for use only inside the ODE solver.

        Parameters
        ----------
        psi : float
            The :math:`\psi` coordinate.

        Returns
        -------
        tuple
            2-tuple containing the calculated derivatives.
        """

    @abstractmethod
    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r"""Returns the electric field strength.

        Input and output are np.ndarrays

        Parameters
        ----------
        psi : np.ndarray
            The ψ values.

        Returns
        -------
        np.ndarray
            1D numpy array with calculated :math:`E` values.

        """

    @abstractmethod
    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r"""Calculates Electric Potential in [V] from ψ.

        Input and output are np.ndarrays

        Parameters
        ----------
        psi : np.ndarray
            The ψ values.

        Returns
        -------

        np.ndarray
            1D numpy array with calculated :math:`\Phi` values.

        """


# =======================================================


class Nofield(ElectricField):
    r"""Initializes an electric field of 0

    Exists to avoid compatibility issues.
    """

    def __init__(self):
        pass

    def e_values(self, psi: float) -> tuple[float, float]:
        return (0, 0)

    def e_valuesNU(self, psi: float) -> tuple[float, float]:
        return (0, 0)

    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        return 0 * psi

    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        return 0 * psi

    def __repr__(self):
        return "No electric field,"


class Radial(ElectricField):
    r"""Initializes an electric field of the form:
    :math:`E(r) = -E_a\exp\bigg[-\dfrac{(r-r_a)^2}{r_w^2})\bigg]`
    """

    def __init__(
        self,
        a: Quantity,
        Ea: Quantity,
        B0: Quantity,
        peak: float,
        rw: float,
    ):
        r"""Initializes the field's parameters.

        Parameters
        ----------

        a : Quantity
            The tokamak's minor radius in [m].
        Ea : Quantity
            The Electric field magnitude in [V/m].
        B0 : Quantity
            The Magnetic field's strength in [T].
        peak : float
            The Electric field's peak point with respect to
            :math:`\psi_{wall}`.
        rw : float
            The Electric field's waist width relative to :math:`r_{wall}`,
            defined as:
            :math:`r_{waist} = \alpha\cdot \text{rw}`.

        """

        # SI Quantities
        self.a = a.to("meters")
        self.Ea = Ea.to("Volts/meter")
        self.B0 = B0.to("Tesla")
        self.rpeak = (self.a * peak).to("meters")
        self.rw = (self.a * rw).to("meters")  # waist, not wall
        self.psia = (self.B0 * self.rpeak**2 / 2).to("Magnetic_flux")
        self.psiw = (self.B0 * self.rw**2 / 2).to("Magnetic_flux")

        # [NU] Conversions for the derivatives:
        self.EaNU = self.Ea.to("NUVolts_per_NUmeter")
        self.psiaNU = self.psia.to("NUMagnetic_flux")
        self.psiwNU = self.psiw.to("NUMagnetic_flux")

        # Unitless quantities, makes it a bit faster if defined here
        self._Ea = self.Ea.magnitude
        self._rpeak = self.rpeak.magnitude
        self._rw = self.rw.magnitude
        self._psia = self.psia.magnitude
        self._psiw = self.psiw.magnitude
        self._sr_psia = sqrt(self._psia)
        self._sr_psiw = sqrt(self._psiw)
        self._EaNU = self.EaNU.magnitude
        self._psiaNU = self.psiaNU.magnitude
        self._sr_psiaNU = sqrt(self._psiaNU)
        self._psiwNU = self.psiwNU.magnitude
        self._sr_psiwNU = sqrt(self._psiwNU)

    def solverE(self, psi: float) -> tuple[float, float]:
        Phi_der_psi = (
            self._EaNU
            / (sqrt(2 * psi))
            * exp(-((sqrt(psi) - self._sr_psiaNU) ** 2) / self._psiwNU)
        )
        Phi_der_theta = 0

        return [Phi_der_psi, Phi_der_theta]

    def Er_of_psi(self, psi: np.ndarray) -> np.ndarray:
        r = np.sqrt(2 * psi)
        E = -self._Ea * np.exp(-((r - self._rpeak) ** 2) / self._rw**2)
        return E

    def Phi_of_psi(self, psi: np.ndarray) -> np.ndarray:
        Phi = (
            self._EaNU
            * np.sqrt(np.pi * self._psiwNU / 2)
            * (
                erf((np.sqrt(psi) - self._sr_psiaNU) / self._sr_psiwNU)
                + erf(self._sr_psiaNU / self._sr_psiwNU)
            )
        )
        return Phi

    def __repr__(self):
        return (
            "Radial: "
            + f"Ea={self.Ea:.4g~P}, peak={self.rpeak:.4g~P}, rw={self.rw:.4g~P}."
        )
