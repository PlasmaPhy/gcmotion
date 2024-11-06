r"""
About Magnetic field objects
----------------------------

A Magnetic Field object is a class instance containing all the information 
about the magnetic field of the system. It implements all the methods needed 
buy the solver and other calculations, and is called automatically wherever required.

To add a new Magnetic Field, simply copy-paste an already existing class
(idealy the Nofield one) and fill the ``__init__()``, ``B()`` and ``B_der()``
methods to fit your Magnetic Field.  In case your Magnetic field has extra 
parameters you want to pass as arguments, you must also create an 
``__init__()`` method and declare them. To avoid errors, your class should
inherit the ``MagneticField`` class.

The general structure is this::

    class MyMagneticField(MagneticField):

        def __init__(self, *parameters):
            "Parameter setup."
    
        def b(self, r, theta): 
            return b

        def b_der(self, r, theta)
            return (B_der_psi, B_der_theta)

.. note::
    The Magnetic Field's attributes must be Quantities with units in SI so 
    they can be referenced from anywhere in the code, but the methods input
    **must** be purely numeric. As a result, the output is also purely 
    numeric and in the same units as the input.
  
.. rubric:: The 'MagneticField' Abstract Base Class

The base class that every other class inherits from is ``MagneticField``. 
This class does nothing, it is only a template.

.. autoclass:: MagneticField
    :member-order: bysource
    :members: __init__, b, b_der

"""

import pint
import numpy as np
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
    def b_values(self, psi: float, theta: float) -> tuple[float, float]:
        r"""Derivatives of :math:`b(\psi, \theta)` with respect to
        :math:`r, \theta`.

        Intended for use only inside the ODE solver.

        Parameters
        ----------
        psi : float
            The :math:`\psi` coordinate.
        theta : float
            The :math:`\theta` coordinate.

        Returns
        -------
        tuple
            2-tuple containing the calculated derivatives.

        """

    @abstractmethod
    def b(
        self, r: float | np.ndarray, theta: float | np.ndarray
    ) -> float | np.ndarray:
        r"""Returns the magnetic field strength.Output size is the same
        as input size.

        When used inside the solver it returns :math:`b`, the strength
        of the magnetic field normalized to its value on magnetic axis.
        When used for the tokamak profile plotting, it returns an np.ndarray.

        Parameters
        ----------
        r : float | np.ndarray
            the r position of the particle.
        theta : float | np.ndarray
            the theta position of the particle

        Returns
        -------
        float | np.ndarray
            The magnetic field strength.

        """


# =======================================================


class LAR(MagneticField):
    r"""Initializes the standard Large Aspect Ration magnetic field."""

    def __init__(self, B0: Quantity, i: Quantity, g: Quantity):
        r"""Parameters initialization.

        Parameters
        -----------
        B0 : Quantity
            The strength of the magnetic field *on the magnetic axis*
            in [T].
        i : Quantity
            The toroidal current in [T * m].
        g : Quantity
            The poloidal current in [T * m].
        """

        # SI Quantities
        self.B0 = B0.to("Tesla")
        self.i = i.to("Plasma_Current")
        self.g = g.to("Plasma_Current")

        # [NU] Conversions
        self.B0NU = B0.to("NUTesla")
        self.iNU = i.to("NUPlasma_Current")
        self.gNU = g.to("NUPlasma_Current")

    def b_values(self, psi: float, theta: float):
        root = sqrt(2 * psi)
        cos_theta = cos(theta)
        b = 1 - root * cos_theta
        b_der_psi = cos_theta / root
        b_der_theta = root * sin(theta)

        return (b, b_der_psi, b_der_theta)

    def b(self, r: float | np.ndarray, theta: float | np.ndarray):
        if isinstance(r, (int, float)):
            return 1 - r * cos(theta)
        else:
            return 1 - r * np.cos(theta)

    def __repr__(self):
        return (
            "LAR: "
            + f"B0={self.B0:.4g~P}, I={self.i:.4g~P}, g={self.g:.4g~P}."
        )
