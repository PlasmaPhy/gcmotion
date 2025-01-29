r"""
Defines the MagneticFieldBase class and all available bfield
configurations.
"""

import os
import pint
import numpy as np
import xarray as xr
from termcolor import colored
from scipy.interpolate import UnivariateSpline, RectBivariateSpline

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
    def bigNU(self, phi: float | np.ndarray, theta: float | np.ndarray) -> float | np.ndarray:
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
    def solverbNU(self, psi: float, theta: float) -> tuple[float, float, float]:
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
            colored("LAR", "light_blue") + f": B0={self.B0:.4g~}, I={self.i:.4g~}, g={self.g:.4g~}."
        )


class Smart(MagneticField):

    def __init__(self, triangularity: str = "positive"):
        r""" """

        if triangularity == "positive":
            path = r"C:\Users\georg\OneDrive\Desktop\My Files\Plasma & Fusion\smart_equil_for_internal_use\normalized_Equil-PT-S2-000021-B-updated-COCOS2.nc"
        elif triangularity == "negative":
            path = r"C:\Users\georg\OneDrive\Desktop\My Files\Plasma & Fusion\smart_equil_for_internal_use\normalized_SMART_NT_Jesus_delta_scans_d_11.nc"
        else:
            print('Triangularity must be either "positive" or "negative"')
            return

        try:
            dataset = xr.open_dataset(path)
            self.dataset = dataset
        except FileNotFoundError:
            raise FileNotFoundError(f"No file found at '{path}'")

        # Extract the arrays
        psi_values = dataset.psi.data
        theta_values = dataset.boozer_theta.data
        b_values = dataset.b_field_norm.data.T
        db_dpsi_values = dataset.db_dpsi_norm.data.T
        db_dtheta_values = dataset.db_dtheta_norm.data.T
        i_values = dataset.I_norm.data
        g_values = dataset.g_norm.data

        # Create splines
        self.b_spline = RectBivariateSpline(
            x=theta_values,
            y=psi_values,
            z=b_values,
        )
        self.db_dpsi_spline = RectBivariateSpline(
            x=theta_values,
            y=psi_values,
            z=db_dpsi_values,
        )
        self.db_dtheta_spline = RectBivariateSpline(
            x=theta_values,
            y=psi_values,
            z=db_dtheta_values,
        )

        self.i_spline = UnivariateSpline(
            x=psi_values,
            y=i_values,
        )
        self.g_spline = UnivariateSpline(
            x=psi_values,
            y=g_values,
        )

        self.ider_spline = self.i_spline.derivative(n=1)
        self.gder_spline = self.g_spline.derivative(n=1)

        # Useful attributes
        _B0 = float(dataset.Baxis.data)  # Tesla
        self.B0 = pint.UnitRegistry.Quantity(_B0, "Tesla")
        self.is_numerical = True

    def bigNU(self, psi: float | np.ndarray, theta: float | np.ndarray):

        theta = theta % (2 * np.pi)
        b = self.b_spline(x=theta, y=psi, grid=False)
        i = self.i_spline(x=psi)
        g = self.g_spline(x=psi)

        return (b, i, g)

    def solverbNU(self, psi: float, theta: float):

        theta = theta % (2 * np.pi)
        # Field and currents
        b, i, g = self.bigNU(psi, theta)

        # Field derivatives
        db_dpsi = self.db_dpsi_spline(x=theta, y=psi, grid=False)
        db_dtheta = self.db_dtheta_spline(x=theta, y=psi, grid=False)

        # Current derivatives
        i_der = self.ider_spline(psi)
        g_der = self.gder_spline(psi)

        # Pack them up
        currents = (i, g)
        b_der = (db_dpsi, db_dtheta)
        currents_der = (i_der, g_der)

        return b, b_der, currents, currents_der

    def __repr__(self):
        return colored("Smart", "light_blue") + f": B0={self.B0:.4g~}."
