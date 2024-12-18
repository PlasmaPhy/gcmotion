import os
import xarray as xr

from gcmotion.utils.quantity_constructor import QuantityConstructor


class SmartInit:
    r"""Imports the necessary constants (R, B0, :math:`\psi_{wall}`) from the
    dataset to initialize the Quantity Constructor and Tokamak

    The particle's species must be defined by the user, while the rest of the
    Quantities are already stored in the dataset.

    Parameters
    ----------
    species : {'p', 'e', 'D', 'T', 'He3', 'He4'}
        The particle's species. This field is case-insensitive.

    """

    def __init__(self, species):

        # Open the dataset
        parent = os.path.dirname(__file__)
        path = os.path.join(parent, "smart.nc")
        try:
            ds = xr.open_dataset(path)
        except FileNotFoundError:
            raise FileNotFoundError(f"No file found at '{path}'")

        # Get numerical values from dataset
        B0 = float(ds.Baxis.data)  # Tesla
        R = ds.raxis.data  # meters
        psi_wallNU = float(ds.psi[-1].data)  # NUMagnetic_flux

        # Create Quantity Constructor
        self.Q = QuantityConstructor(
            R=R,
            B0=B0,
            _psi_wallNU=psi_wallNU,
            species=species,
        )

        # Create R, B0 and (phony) a Quantities
        self.R = self.Q(R, "meters")
        self.B0 = self.Q(B0, "Tesla")
        self.psi_wallNU = self.Q(psi_wallNU, "NUMagnetic_flux")
        # self.a = self.Q(
        #     ds.R.sel(psi=1000000, boozer_theta=0, method="nearest").data - R,
        #     "meters",
        # )
        psi_wall = self.psi_wallNU.to("Tesla * meters^2")
        self.a = (2 * psi_wall / self.B0) ** (1 / 2)

    def get_QuantityConstructor(self):
        return self.Q
