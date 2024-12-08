import numpy as np
import matplotlib.pyplot as plt
from typing import Literal

from gcmotion.entities.profile import Profile
from gcmotion.plot._base._base_profile_contour import _base_profile_contourE
from gcmotion.configuration.scripts_configuration import (
    FrequencyConfig,
    QKineticConfig,
)

type flux_units = Literal["Magnetic_flux", "NUMagnetic_flux", "psi_wall"]
type E_units = Literal["eV", "keV", "NUJoule"]


class FrequencyAnalysis:

    def __init__(self, profile: Profile):

        # Grab configuration
        config = FrequencyConfig()
        self.profile = profile

        # Create a phony contour object to calculate a base energy span between
        # psilim = [0, 1.2] (psiwall)
        phony_fig, phony_ax = plt.subplots()
        autoC = _base_profile_contourE(
            profile=self.profile,
            ax=phony_ax,
            thetalim=config.auto_thetalim,
            psilim=config.auto_pslim,
            levels=config.auto_levels,
            mode="lines",
        )
        del phony_fig, phony_ax
        plt.close()

        self.Emin = autoC.zmin
        self.Emax = autoC.zmax
        self._psilim = (autoC._mins[1], autoC._maxs[1])

        pass

    def theta_frequencies(self, **args):
        r"""

        Parameters
        ----------
        psilim : list, optional
            The psi span in which to calculate :math:`\hat\omega_\theta`.
            Defaults to [1, 1.2].
        levels : int, optional
            The Energy levels to calculate the corresponding frequency.
            Defaults to 250.
        grid_density: int, optional
            The contouring grid density. Defaults to 200.

        Other Parameters
        ----------------
        flux_units : {"Magnetic_flux", "NUMagnetic_flux", "psi_wall"}, optional
            The :math:`\psi, P_\theta` units of the contour. Defaults to
            "Magnetic_flux" ("Tesla * meters^2")
        E_units : {"eV", "keV", "NUJoule"}, optional
            The Energy units. Defaults to "keV".
        potential : bool, optional
            Whether or not to add the Electric Potential :math:`\Phi` term in
            the Hamiltonian. Defaults to True.
        wall : bool, optional
            Whether or not to shade the area :math:`\psi > \psi_{wall}`.
            Defaults to True.
        """

        pass

    def zeta_frequencies():
        pass

    def qkinetic(self, **args):
        r"""

        Parameters
        ----------
        Espan : list of Quantities
            The Energy span to calculate qkinetic.

        """

        # Grab configuration and overwrite passed values
        config = QKineticConfig()
        for key, value in args.items():
            setattr(config, key, value)

        Emin = args.get("Emin", self.Emin)
        Emax = args.get("Emax", self.Emax)
        Espan = np.linspace(Emin, Emax, config.levels)

        omega_thetas = self.omega_thetas(Espan=Espan)
        omega_zetas = self.omega_zetas(Espan=Espan)

        pass

    def freq_from_area(self):
        pass
