r"""
Defines the MagneticFieldBase class and all available bfield
configurations.
"""

import os
import pint
import numpy as np
import xarray as xr
from termcolor import colored
from scipy.interpolate import RectBivariateSpline, make_interp_spline
from time import time

from gcmotion.utils.logger_setup import logger
from gcmotion.configuration.scripts_configuration import (
    NumericalDatasetsConfig,
)

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
    def solverbNU(self, psi: float, theta: float) -> tuple[float]:
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


class NumericalMagneticField(MagneticField):
    r"""Numerical Magnetic field base class

    Opens the dataset and creates the splines needed for the querry methods.
    Also defines Bmin and Bmax, used in the parabolas.

    ``solverB`` and ``bigNU`` work identically to the analytic Magnetic field
    methods.

    Parameters
    ----------
    filename : str
        The "\*.nc" file located at gcmotion/tokamak/reconstructed.

    Attributes
    ----------
    Bmin, Bmax : Quantities
        The minimum and maximum magnetic field values.
    theta_min, theta_max : Quantities
        The :math:`\theta` coordinates where the magnetic field takes its
        minimum/maximum value.
    psi_min, psi_max : Quantities
        The :math:`\psi` coordinates where the magnetic field takes its
        minimum/maximum value.

    """

    is_analytical = False
    is_numerical = True

    def __init__(self, filename: str):

        # Open the dataset
        parent = os.path.dirname(__file__)  # Relative to this directory
        path = os.path.join(parent, "reconstructed", filename)
        try:
            dataset = xr.open_dataset(path)
            self.dataset = dataset
            logger.info("Dataset initialized correctly.")
        except FileNotFoundError:
            logger.error("Error opening Dataset.")
            raise FileNotFoundError(f"No file found at '{path}'")

        # Extract the arrays
        interval = NumericalDatasetsConfig.boozer_theta_downsampling_factor
        psi_values = dataset.psi.data
        theta_values = dataset.boozer_theta.data[::interval]
        b_values = dataset.b_field_norm.data.T[::interval]
        i_values = dataset.I_norm.data
        g_values = dataset.g_norm.data

        # Extrapolate psi to containn psi=0
        psi_values = np.insert(psi_values, 0, 0)

        # Extrapolate all arrays so their value at the axis is the same as
        # their value on the closest point to the axis to avoid errors when
        # integrating.
        i_values = np.insert(i_values, 0, i_values[0])
        g_values = np.insert(g_values, 0, g_values[0])
        axis_values = np.full((b_values.shape[0], 1), 1)
        b_values = np.hstack((axis_values, b_values))

        # Create splines
        self.b_spline = RectBivariateSpline(
            x=theta_values,
            y=psi_values,
            z=b_values,
        )
        # NOTE: Do not use the dataset bfield derivatives. They introduce small
        # non-Hamiltonian terms, resulting in a noticeable fluxuation and even
        # loss of Energy, due to error propagation through 2 different splines.
        # Use the derivatives calculated from the bfield spline instead.
        self.db_dpsi_spline = self.b_spline.partial_derivative(
            dx=0,
            dy=1,
        )
        self.db_dtheta_spline = self.b_spline.partial_derivative(
            dx=1,
            dy=0,
        )
        self.i_spline = make_interp_spline(
            x=psi_values,
            y=i_values,
            k=NumericalDatasetsConfig.currents_spline_order,
        )
        self.g_spline = make_interp_spline(
            x=psi_values,
            y=g_values,
            k=NumericalDatasetsConfig.currents_spline_order,
        )

        self.ider_spline = self.i_spline.derivative(nu=1)
        self.gder_spline = self.g_spline.derivative(nu=1)

        # Useful attributes
        Q = pint.get_application_registry().Quantity
        _B0 = float(dataset.Baxis.data)  # Tesla
        self.B0 = Q(_B0, "Tesla")

        # Magnetic field strength extremum
        start = time()
        da = dataset.b_field_norm
        mins = da.where(da == da.min(), drop=True).squeeze()
        maxs = da.where(da == da.max(), drop=True).squeeze()
        # Flatten()[0] is needed since it sometimes finds 2 extremums with
        # effectively the same value
        self.Bmin = Q(mins.values.flatten()[0], "NUTesla").to("Tesla")
        self.Bmax = Q(maxs.values.flatten()[0], "NUTesla").to("Tesla")
        self.theta_min = float(mins.boozer_theta.values.flatten()[0])
        self.theta_max = float(maxs.boozer_theta.values.flatten()[0])
        self.psi_min = Q(mins.psi.values, "NUMagnetic_flux")
        self.psi_max = Q(maxs.psi.values, "NUMagnetic_flux")
        end = time()
        duration = Q(end - start, "seconds")
        logger.info(f"Numerical bfield extremum search took {duration:.4g~#P}")

    def bigNU(
        self, psi: float | np.ndarray, theta: float | np.ndarray
    ) -> tuple:
        theta = theta % (2 * np.pi)
        b = self.b_spline(x=theta, y=psi, grid=False)
        i = self.i_spline(x=psi)
        g = self.g_spline(x=psi)

        if isinstance(psi, (float, int)):
            return (float(b), float(i), float(g))
        else:
            return (b, i, g)

    def solverbNU(self, psi: float, theta: float) -> tuple:

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
        b = float(b)
        currents = (float(i), float(g))
        b_der = (float(db_dpsi), float(db_dtheta))
        currents_der = (float(i_der), float(g_der))

        return b, b_der, currents, currents_der


# ============================================================================


class LAR(MagneticField):
    r"""Initializes the standard Large Aspect Ratio magnetic field.

    Parameters
    ----------
    B0 : Quantity
        The strength of the magnetic field *on the magnetic axis*
        in [T].
    i : Quantity
        The toroidal current in units of [Plasma Current].
    g : Quantity
        The poloidal current in units of [Plasma Current].

    """

    is_analytical = True
    is_numerical = False

    def __init__(self, B0: Quantity, i: Quantity, g: Quantity):

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
        self.plain_name = "LAR"

        # Minimum/Maximum values and locations
        Q = pint.get_application_registry().Quantity
        self.psi_wall = Q(1, "psi_wall")
        self.theta_min, self.theta_max = 0, np.pi
        self.psi_min = self.psi_max = self.psi_wall.to("NUMagnetic_flux").m

        _BminNU = self.bigNU(psi=self.psi_min, theta=self.theta_min)[0]
        _BmaxNU = self.bigNU(psi=self.psi_max, theta=self.theta_max)[0]
        self.Bmin = Q(_BminNU, "NUTesla").to("Tesla")
        self.Bmax = Q(_BmaxNU, "NUTesla").to("Tesla")

    def bigNU(
        self, psi: float | np.ndarray, theta: float | np.ndarray
    ) -> tuple[float]:

        if isinstance(psi, float):
            b = 1 - sqrt(2 * psi) * cos(theta)
            g = self._gNU
            i = self._iNU
        elif isinstance(psi, np.ndarray):
            b = 1 - np.sqrt(2 * psi) * np.cos(theta)
            g = np.tile(self._gNU, psi.shape)
            i = np.tile(self._iNU, psi.shape)

        return (b, i, g)

    def solverbNU(self, psi: float, theta: float) -> tuple:

        # Field and currents
        b, i, g = self.bigNU(psi, theta)

        # Field derivatives
        root = sqrt(2 * psi)
        cos_theta = cos(theta)
        b_der_psi = -cos_theta / root
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


class SmartPositive(NumericalMagneticField):
    r"""Initializes a bfield object with numerical data from the Smart Tokamak
    with **Positive** Triangularity.

    The dataset must be stored in
    *./gcmotion/tokamak/reconstructed/smart_positive.nc*.

    """

    def __init__(self):
        filename = "smart_positive.nc"
        super().__init__(filename=filename)

        self.plain_name = "Smart - Positive"

    def __repr__(self):
        return (
            colored("Smart - Positive", "light_blue") + f": B0={self.B0:.4g~}."
        )


class SmartNegative(NumericalMagneticField):
    r"""Initializes a bfield object with numerical data from the Smart Tokamak
    with **Negative** Triangularity.

    The dataset must be stored in
    *./gcmotion/tokamak/reconstructed/smart_negative.nc*.

    """

    def __init__(self):
        filename = "smart_negative.nc"
        super().__init__(filename=filename)

        self.plain_name = "Smart - Negative"

    def __repr__(self):
        return (
            colored("Smart - Negative", "light_blue") + f": B0={self.B0:.4g~}."
        )


class SmartNegative2(NumericalMagneticField):
    r"""Initializes a bfield object with numerical data from the Smart Tokamak
    with **Negative** Triangularity.

    The dataset must be stored in
    *./gcmotion/tokamak/reconstructed/smart_negative.nc*.

    """

    def __init__(self):
        filename = "smart_negative2.nc"
        super().__init__(filename=filename)

        self.plain_name = "Smart - Negative 2"

    def __repr__(self):
        return (
            colored("Smart - Negative 2", "light_blue")
            + f": B0={self.B0:.4g~}."
        )


class DTTPositive(NumericalMagneticField):
    r"""Initializes a bfield object with numerical data from the Divertor Test
    Tokamak with **Positive** Triangularity..

    The dataset must be stored in
    *./gcmotion/tokamak/reconstructed/dtt_positive.nc*.

    """

    def __init__(self):
        filename = "dtt_positive.nc"
        super().__init__(filename=filename)

        self.plain_name = "DTT - Positive"

    def __repr__(self):
        return (
            colored("DTT - Positive", "light_blue") + f": B0={self.B0:.4g~}."
        )


class DTTNegative(NumericalMagneticField):
    r"""Initializes a bfield object with numerical data from the Divertor Test
    Tokamak with **Negative** Triangularity.

    The dataset must be stored in
    *./gcmotion/tokamak/reconstructed/dtt_negative.nc*.

    """

    def __init__(self):
        filename = "dtt_negative.nc"
        super().__init__(filename=filename)

        self.plain_name = "DTT - Negative"

    def __repr__(self):
        return (
            colored("DTT - Negative", "light_blue") + f": B0={self.B0:.4g~}."
        )
