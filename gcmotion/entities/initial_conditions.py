from pint import Quantity
from termcolor import colored

from gcmotion.entities.profile import Profile
from gcmotion.configuration.physical_constants import PhysicalConstants

from gcmotion.utils.logger_setup import logger

type QuantityArray = Quantity


class InitialConditions:
    r"""Returns an InitialConditions entity.

    This entity's only atrributes are its given arguements and their
    corresponding values in NU.

    Parameters
    ----------
    theta0, zeta0 : float | Quantity
        The particle's initial :math:`\theta_0` and :math:`\zeta_0`.
    psi0: Quantity
        The :math:`\psi_0` initial values in dimensions of [Magnetic_flux] or
        [psi_wall]. [psi_wall] is a unit of Magnetic flux, where 1[psi_wall] is
        defined to be the Magnetic flux of the last closed surface. This way we
        can set the :math:`\psi_0` initial value with respect to the tokamak's
        :math:`\psi_{wall}`.
    t_eval : Quantity (1D array)
        The time interval return values, [:math:`t_0, t_f, steps`],
        in dimensions of [time].
    muB : Quantity, optional
        The particle's initial :math:`\mu B` value in units of [energy]. If
        specified, this value overwrites the `mu` attribute in the
        PhysicalParameters object.

    Notes
    -----
    An InitialConditions object contains the following attributes, which
    include the input arguements:

        #. theta0, zeta0 : float | Quantity
            The particle's initial :math:`\theta_0` and :math:`\zeta_0`.
        #. psi0, psi0NU : Quantity
            The initial :math:`\psi_0` value in SI/NU.
        #. t_eval, t_evalNU : QuantityArray
            The time interval return values in SI/NU.
        #. muB, muBNU : None | Quantity
            The particle's initial :math:`\mu B` value in SI/NU.

    Example
    -------
    How to initialize an `InitialConditions` object.

    >>> import gcmotion as gcm
    >>> import numpy as np
    >>>
    >>> # Quantity Constructor
    >>> Rnum = 1.65
    >>> anum = 0.5
    >>> B0num = 1
    >>> species = "p"
    >>> ureg, Q = gcm.setup_pint(R=Rnum, a=anum, B0=B0num, species=species)
    >>>
    >>> # Intermediate values
    >>> R = Q(Rnum, "meters")
    >>> a = Q(anum, "meters")
    >>> B0 = Q(B0num, "Tesla")
    >>> i = Q(0, "NUPlasma_current")
    >>> g = Q(1, "NUPlasma_current")
    >>> Ea = Q(73500, "Volts/meter")
    >>>
    >>> init = gcm.InitialConditions(
    ...     muB=Q(0.5, "keV"),
    ...     theta0=0,
    ...     zeta0=0,
    ...     psi0=Q(0.6, "psi_wall"),
    ...     t_eval=Q(np.linspace(0, 1e-3, 1000), "seconds"),
    ... )

    """

    def __init__(
        self,
        theta0: float | Quantity,
        zeta0: float | Quantity,
        psi0: Quantity,
        t_eval: QuantityArray,
        muB: Quantity = None,
    ):
        r"""Initializes parameters and corresponding NU Quantities"""
        logger.info("==> Initializing InitialConditions...")

        # muB parsing
        # The actual decision is made in the Particle class.
        self.muB = None if muB is None else muB.to("keV")

        # Keep the angles' type
        self.theta0 = float(theta0)
        self.zeta0 = float(zeta0)

        self.psi0 = psi0.to("Magnetic_flux")
        self.t_eval = t_eval.to("seconds")

        # Corresponding NU values
        self.psi0NU = self.psi0.to("NUMagnetic_flux")
        self.t_evalNU = self.t_eval.to("NUseconds")

        logger.info("\t InitialConditions Initialization Complete")

    def __repr__(self):
        t0 = self.t_eval[0].m
        dt = self.t_eval[1].m - t0
        tf = self.t_eval[-1].m
        muB = f"{self.muB:.4g~P}" if isinstance(self.muB, Quantity) else None

        return (
            f"theta0 = {self.theta0:.4g}, "
            + f"zeta0 = {self.zeta0:.4g}, "
            + f"psi0 = {self.psi0:.4g~P}, "
            + f"t_eval: [t0, tf, dt]=[{t0:.4g}, {tf:.4g}, {dt:.4g}]s, "
            + f"muB = {muB}, "
        )

    def __str__(self):
        t0 = self.t_eval[0].m
        dt = self.t_eval[1].m - t0
        tf = self.t_eval[-1].m
        muB = f"{self.muB:.4g~P}" if isinstance(self.muB, Quantity) else None

        return (
            colored("\nInitial Conditions :\n", "green")
            + f"{"theta0":>23} : {f'{self.theta0:.4g} radians':<16} \n"
            + f"{"zeta0":>23} : {f'{self.zeta0:.4g} radians':<16} \n"
            + f"{"psi0":>23} : {f'{self.psi0:.4g~}':<16}"
            + f"({self.psi0NU:.4g})\n"
            + f"{"t_eval":>23} : "
            + f"{f'[t0, tf, dt] = [{t0:.4g}, {tf:.4g}, {dt:.4g}]s':<16}\n"
            + f"{"muB":>23} : {muB}"
        )


class _InitialConditionsFull(InitialConditions, Profile):
    r"""Extends the InitialConditions class to include the extra initial
    conditions :math:`\rho_0`, :math:`\psi_{p0}` and :math:`P_\theta_0`
    """

    def __init__(self, init, profile):
        r"""Grabs attributes from parent classes and calculates the extra
        initial conditions."""
        logger.debug("==> Initializing InitialConditionFull...")

        # Grab attributes from parents
        InitialConditions.__init__(
            self,
            theta0=init.theta0,
            zeta0=init.zeta0,
            psi0=init.psi0,
            t_eval=init.t_eval,
            muB=init.muB,
        )

        Profile.__init__(
            self,
            tokamak=profile.tokamak,
            params=profile.params,
        )

        # Grab particle's mass and charge
        M = getattr(PhysicalConstants, profile.species + "_M")
        Z = getattr(PhysicalConstants, profile.species + "_Z")

        self.miNU = self.Q(M, "Proton_mass")
        self.qiNU = self.Q(Z, "Proton_charge")

        self.mi = self.miNU.to("kilogram")
        self.qi = self.qiNU.to("Coulomb")
        # Initial (B,i,g) and Phi values
        # Note: bigNU() and PhiNU() input and output are in [NU]
        # Symbolize the magnitude of the repspective variable with a _
        _psi0NU = self.psi0NU.m
        _theta0 = self.theta0
        _b0NU, _i0NU, _g0NU = self.bfield.bigNU(_psi0NU, _theta0)
        _Phi0NU = float(self.efield.PhiNU(_psi0NU, _theta0))
        B_init = (self.B0 * _b0NU).to("Tesla")
        i_init = self.Q(_i0NU, "NUPlasma_current").to("Plasma_current")
        g_init = self.Q(_g0NU, "NUPlasma_current").to("Plasma_current")
        Phi_init = self.Q(_Phi0NU, "Volts")

        # muB parsing
        # If "muB" is not present in "init", use "params.mu" to set mu and muB.
        # Otherwise, use "init.muB",
        # recalculate mu and update particle's profile.
        if init.muB is None:
            self.mu = profile.mu.to("kev/T")
            self.muB = (profile.mu * B_init).to("keV")
            mu_selection = "profile"
        else:
            # parse muB depending on its dimensionality
            if (
                init.muB.dimensionality
                == self.Q("Magnetic_moment").dimensionality
            ):
                self.mu = init.muB.to("keV/Tesla")
                self.muB = (self.mu * B_init).to("keV")
            elif init.muB.dimensionality == self.Q("keV").dimensionality:
                self.muB = init.muB.to("keV")
                self.mu = (self.muB / B_init).to("keV/Tesla")
            mu_selection = "init"
        logger.debug(f"\tGrabbed mu from {mu_selection}.")
        logger.debug(f"\tmu={self.mu:.4g~P}, muB={self.muB:.4g~P}")

        # Derivative Quantities
        self.psip0 = self.Q(  # psipNU only accepts [NU]
            self.qfactor.psipNU(_psi0NU),
            "NUMagnetic_flux",
        ).to("Magnetic_flux")

        self.rho0 = ((self.Pzeta + self.psip0) / g_init).to(
            "meters"
        )  # Note: [rho] = Magnetic_flux / Plasma_current = meters)

        self.Ptheta0 = (self.psi0 + self.rho0 * i_init).to("Magnetic_flux")

        # Energies
        self.E = (
            self.qi**2 / (2 * self.mi) * self.rho0**2 * B_init**2
            + self.mu * B_init
            + self.qi * Phi_init
        ).to("Joule")

        self.EkeV = self.E.to("keV")

        # Corresponding NU attributes
        self.muNU = self.mu.to("NUMagnetic_moment")
        self.muBNU = self.muB.to("NUJoule")
        self.psip0NU = self.psip0.to("NUMagnetic_flux")
        self.rho0NU = self.rho0.to("NUMeters")
        self.Ptheta0NU = self.Ptheta0.to("NUMagnetic_flux")
        self.ENU = self.E.to("NUJoule")

        logger.debug("\t _InitialConditionsFull initialization completed")
