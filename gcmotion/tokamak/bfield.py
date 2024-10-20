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

.. caution::
    **All values, both input and output are in specific units.**

    Specifically, g and I in [NU] and B in [T].
  
The general structure is this::

    class MyMagneticField(MagneticField):

        def __init__(self, *<parameters>):
            self.id = "foo" # Simple id used only for logging.
            self.params = {} # Tweakable parameters, used only for logging.
            <set parameters>
    
        def B(self, r, theta): 
            if type(r) is float:
                return 1 - r * cos(theta)
            else:
                return 1 - r * np.cos(theta)

        def B_der(self, r, theta)
            return (B_der_psi, B_der_theta)

.. note::
    These two methods return the same type as the input (either Python floats 
    or np.ndarrays). When used inside the solver, they should return a 
    Python float, and not a np.float. Solvers need to be fast so they work 
    with built-in floats, while plotting functions work with np.ndarrays.This 
    is mainly for optimization reasons and should probably not cause problems.
  
.. rubric:: The 'MagneticField' Abstract Base Class

The base class that every other class inherits from is ``MagneticField``. 
This class does nothing, it is only a template.

.. autoclass:: MagneticField
    :members: __init__, B, B_der

"""

import numpy as np
from math import cos, sin, sqrt
from abc import ABC, abstractmethod


class MagneticField(ABC):
    r"""Electric field base class"""

    def __init__(self):
        r"""Contains all the needed parameters"""
        self.id = "Base Class"
        self.params = {}

    @abstractmethod
    def B(self, r: float | np.ndarray, theta: float | np.ndarray) -> float | np.ndarray:
        r"""Returns the magnetic field strength.

        Used inside the solver (input and return type must be float),
        and the countour plot of the magnetic field profile (input and
        return type must be np.array).

        Parameters
        ----------

        r : float | np.ndarray
            the r position of the particle.
        theta : float | np.ndarray
            the theta position of the particle

        Returns
        -------

        B : float | np.ndarray
            The magnetic field strength.

        """
        pass

    @abstractmethod
    def B_der(self, psi: float, theta: float) -> tuple[float, float]:
        r"""Derivaties of B(:math:`\psi`, :math:`\theta`) with respect to r,
        :math:`\theta`.

        Intended for use only inside the ODE solver. Returns the derivatives
        with B in normalized units, so no conversion is needed.


        Parameters
        ----------

        psi : float
            The :math:`\psi` coordinate
        theta : float
            The :math:`\theta` coordinate.

        Returns
        -------

        tuple : 2-tuple of floats
            2-tuple containing the calculated derivatives.

        """


# =======================================================


class LAR(MagneticField):
    r"""Initializes the standard Large Aspect Ration magnetic field."""

    def __init__(self, i: float, g: float, B0: float):
        r"""Parameters initialization.

        Parameters
        -----------

        i : float
            The toroidal current in [NU].
        g : float
            The poloidal current in [NU].
        B0 : float
            The magnetic field strength on the magnetic axis in [T].

        """
        self.id = "LAR"
        self.params = {"I": i, "g": g, "B0": B0}

        self.I, self.g, self.B0 = i, g, B0
        self.is_lar = True

    def B(self, r: float | np.ndarray, theta: float | np.ndarray):
        if isinstance(r, (int, float)):
            return 1 - r * cos(theta)
        else:
            return 1 - r * np.cos(theta)

    def B_der(self, psi: float, theta: float):
        root = sqrt(2 * psi)
        B_der_psi = cos(theta) / root
        B_der_theta = root * sin(theta)

        return (B_der_psi, B_der_theta)
