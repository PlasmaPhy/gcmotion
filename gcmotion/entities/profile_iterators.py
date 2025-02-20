import numpy as np

from pint import Quantity
from numpy.typing import ArrayLike
from gcmotion.entities.profile import Profile


class IterProfilePzeta(Profile):
    r"""Iterator that returns a new profile with updated Pzeta.

    Parameters
    ----------
    pzetaspan: Quantity(ArrayLike) | ArrayLike
        If a unitless array is passed, then it is interpreted as values
        relative to psip_wall. If Quantity(ArrayLike) is passed, then the
        values are used as is.

    """

    def __init__(
        self, pzetaspan: Quantity(ArrayLike) | ArrayLike, *args, **kwargs
    ):
        r"""Instantiates parent Profile and sets up pzeta iterator."""

        # Set the initial Pzeta to None for consistency
        super().__init__(*args, **kwargs | {"Pzeta": None})

        # Define pzetaspan values in NU
        if isinstance(pzetaspan, Quantity):
            self.pzetaspan = pzetaspan.to("NUcanonical_Momentum").m
        else:
            psip_wall = self.tokamak.psip_wallNU.magnitude
            self.pzetaspan = psip_wall * np.array(pzetaspan)

        self.iter_pzeta = np.nditer(self.pzetaspan)

    def __iter__(self):
        return self

    def __next__(self):
        self.PzetaNU = self.Q(next(self.iter_pzeta), "NUCanmom")
        return self


class IterProfileMu(Profile):
    r"""Iterator that returns a new profile with updated mu.

    Parameters
    ----------
    muspan : Quantity(ArrayLike)
        The muspan values, as a Quantity array.

    """

    def __init__(self, muspan: Quantity(ArrayLike), *args, **kwargs):
        r"""Instantiates parent Profile and sets up pzeta iterator."""

        # Set the initial mu to None for consistency
        super().__init__(*args, **kwargs | {"mu": None})
        self.iter_mu = np.nditer(muspan)

    def __iter__(self):
        return self

    def __next__(self):
        self.muNU = self.Q(next(self.iter_mu), "NUMagnetic_moment")
        return self


class IterProfileEnergy(Profile):
    r"""Iterator that returns a new profile with updated Energy.

    Parameters
    ----------
    energyspan : Quantity(ArrayLike)
        The energyspan values, as a Quantity array.

    """

    def __init__(self, energyspan: Quantity(ArrayLike), *args, **kwargs):
        r"""Instantiates parent Profile and sets up pzeta iterator."""

        # Set the initial E to None for consistency
        super().__init__(*args, **kwargs | {"E": None})
        self.iter_energy = np.nditer(energyspan)

    def __iter__(self):
        return self

    def __next__(self):
        self.ENU = self.Q(next(self.iter_energy), "NUJoule")
        return self
