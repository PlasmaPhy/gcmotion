import numpy as np
from time import time
from math import sqrt

from gcmotion.tokamak.efield import Nofield

from gcmotion.utils._logger_setup import logger

from gcmotion.configuration.physical_constants import physical_constants

from gcmotion.scripts.orbit import orbit

from gcmotion.classes.parabolas import Construct


class Particle:
    r"""
    Instanciates a particle.

    The particles holds its properties such as mass, species, etc as well as all its calculated
    quantities as its attributes. We can view them at any time:

    .. code-block:: python

        >>> particle.species
        'e'
        >>> particle.Ptheta
        array([0.02601457, 0.02580392, 0.02559509, ..., 0.03456111, 0.03481527,
            0.0350704 ])

    Physical properties and tokamak configuration are automatically setup upon the particle's
    initialization, and its orbit is calculated upon calling its
    :py:meth:`~Particle.run` method.

    Here is a list of the most important attributes:

    #. Initial conditions:
        theta0, psi0, zeta0, rho0, psip0, Ptheta0, Pzeta0
    #. Constants of motion:
        E_NU, E_J, E_eV, mu, Pzeta
    #. Time evolution arrays:
        theta, psi, zeta, rho, psip, Ptheta, Pzeta
    #. Configuration objects and parameters:
        R, a, qfactor, Bfield, Efield, psi_wall, psip_wall
    #. Conversion factors and physical properties
        mass_amu, mass_kg, w0, E_unit,
        Volts_to_NU, NU_to_J, NU_to_eV

    .. hint::
        We can always view all the currently stored attributes by running:

        .. code-block:: python

            >>> from pprint import pprint
            >>> pprint(vars(cwp))

    Example
    -------

    See :ref:`Creating a particle <user_guide_particle_creation>`

    .. rubric:: Methods
        :heading-level: 1
    """

    def __init__(
        self,
        tokamak: dict,
        parameters: dict,
    ):
        r"""Initializes particle and grabs configuration.

        :meta public:

        Parameters
        ----------

        tokamak : dict
            A dict containing the tokamak's configuration:

                R : float
                    The tokamak's major radius in [m].
                a : float
                    The tokamak's minor radius in [m].
                qfactor : :py:class:`~gcmotion.tokamak.qfactor.QFactor`
                    Qfactor object that supports query methods for getting values
                    of :math:`q(\psi)` and :math:`\psi_p(\psi)`.
                Bfield : :py:class:`~gcmotion.tokamak.bfield.MagneticField`
                    Magnetic Field Object that supports query methods for getting values of the
                    field magnitude and its derivatives.
                Efield : :py:class:`~gcmotion.tokamak.efield.ElectricField`
                    Electric Field Object that supports query methods for getting values of the
                    field itself and the derivatives of its potential.
        parameters : dict
            A dict containing all the particle-specific parameters.

            t_eval : np.array
                The ODE time interval return values, [:math:`t_0, t_f`, steps], in normalised
                time units [NU].
            theta0, psi0, zeta0, Pzeta0 : floats
                The particle's initial conditions:
                [:math:`\theta_0, \psi_0, \zeta_0, P_{\zeta_0}`]) [NU]
            mu : float
                The magnetic moment in [NU].
            species : str
                The particle species, used to later set charge and mass automatically
                (from :py:mod:`~gcmotion.configuration.physical_constants`)

        """
        logger.info("--------Initializing particle--------")

        def setup_constants():
            """Grabs particle's constants from ``physical_constants.py``"""

            logger.info("Setting up particle's constants...")

            self.species = parameters["species"].lower()
            self.mNU = physical_constants[self.species + "_mNU"]
            self.qNU = physical_constants[self.species + "_qNU"]
            self.mp = physical_constants["mp"]

            logger.debug(f"\tParticle is of species '{self.species}'.")
            logger.info("--> Particle's constants setup successful")

        def setup_tokamak():
            """Sets up tokamak-related attributes."""

            logger.info("Setting up Tokamak...")

            # Dimensions
            self.R, self.a = tokamak["R"], tokamak["a"]

            logger.debug(f"\tTokamak dimensions: R = {self.R}, a = {self.a}")

            # Objects
            self.qfactor = tokamak["qfactor"]
            self.Bfield = tokamak["Bfield"]
            self.Bfield.B0NU = self.mNU / self.qNU
            Efield = tokamak["Efield"]
            if Efield is None or isinstance(Efield, Nofield):
                self.Efield = Nofield()
                self.has_efield = False
            else:
                self.Efield = Efield
                self.has_efield = True

            logger.debug(
                f"\t'{self.qfactor.id}' qfactor used with parameters {self.qfactor.params}"
            )
            logger.debug(
                f"\t'{self.Bfield.id}' Bfield used with parameters {self.Bfield.params}"
            )
            logger.debug(
                f"\t'{self.Efield.id}' Efield used with parameters {self.Efield.params}"
            )

            self.psi_wall = (self.a) ** 2 / 2  # normalized to R
            self.psip_wall = self.qfactor.psip_of_psi(self.psi_wall)

            # psi_p > 0.5 warning
            if self.psip_wall >= 0.5:
                logger.warning(
                    f"\tWARNING: psip_wall = {self.psip_wall:.5g} >= 0,5."
                    + "Parabolas and other stuff will probably not work"
                )

            logger.debug(
                f"\tDerivative quantities: psi_wall = {self.psi_wall:.5g}"
                + f" (CAUTION: normalised to R), psip_wall = {self.psip_wall:.5g}"
            )

            logger.info("--> Tokamak setup successful.")

        def setup_parameters():
            """Sets up the particles initial condition and parameters, as well as the solver's S0."""

            logger.info("Setting up particle's initial conditions...")

            self.t_eval = parameters["t_eval"]
            self.mu = parameters["mu"]
            self.theta0 = parameters["theta0"]
            self.psi0 = (
                parameters["psi0"] * self.psi_wall
            )  # CAUTION! Normalize it to psi_wall
            self.zeta0 = parameters["zeta0"]
            self.Pzeta0 = parameters["Pzeta0"]
            self.psip0 = self.qfactor.psip_of_psi(self.psi0)
            self.rho0 = (
                self.Pzeta0 + self.psip0
            ) / self.Bfield.g  # Pz0 + psip0
            self.Ptheta0 = self.psi0 + self.rho0 * self.Bfield.I  # psi + rho*I

            logger.debug(
                "ODE initial conditions:\n"
                + f"\ttheta0 = {self.theta0:.5g}, psi0 = {self.psi0:.5g}, zeta0 = {self.zeta0:.5g}, Pzeta0 = {self.Pzeta0:.5g}."
            )
            logger.debug(
                "\tOther initial conditions:\n"
                + f"\tPtheta0 = {self.Ptheta0:.5g}, psip0 = {self.psip0:.5g}, rho0 = {self.rho0:.5g}, mu = {self.mu:.5g}"
            )
            logger.debug(
                f"\tTime span (t0, tf, steps): ({self.t_eval[0]}, {self.t_eval[-1]}, {len(self.t_eval)})"
            )
            logger.info("--> Initial conditions setup successful.")

        def setup_logic_flags():
            """Sets up logic flags and initializes variables that must have an initial value"""

            logger.info("Setting up logic flags...")

            self.calculated_conversion_factors = False
            self.calculated_energies = False
            self.calculated_orbit_type = False
            self.calculated_orbit = False
            self.t_or_p = "Unknown"
            self.l_or_c = "Unknown"
            self.percentage_calculated = 0

            # Stored initially to avoid attribute errors
            self.z_0freq = self.z_freq = self.theta_0freq = self.theta_freq = (
                None
            )

            logger.info("--> Logic flags setup successful.")

        setup_constants()
        setup_tokamak()
        setup_parameters()
        setup_logic_flags()

        logger.info("--------Particle Initialization Completed--------\n")

    def __str__(self):
        info_str = (
            "Constants of motion:\n"
            + "\tParticle Energy (normalized):\tE = {:e}\n".format(self.E_NU)
            + "\tParticle Energy (eV):\t\tE = {:e} eV\n".format(self.E_eV)
            + "\tParticle Energy (J):\t\tE = {:e} J\n".format(self.E_J)
            + f"\tToroidal Momenta:\t\tPζ = {self.Pzeta0}\n\n"
            + "Other Quantities:\n"
            + f'\tParticle of Species:\t\t"{self.species}"\n'
            + f"\tOrbit Type:\t\t\t{self.orbit_type_str}\n"
            + f"\tMajor Radius:\t\t\tR = {self.R} meters\n"
            + f"\tMinor Radius:\t\t\tα = {self.a} meters\n"
            + "\tToroidal Flux at wall:\t\tψ = {:n}\n".format(self.psi_wall)
            + "\tTime unit:\t\t\tω = {:e} Hz \n".format(self.w0)
            + "\tEnergy unit:\t\t\tE = {:e} J \n\n".format(self.E_unit)
            + self.solver_output
        )

        return info_str

    def __repr__(self):
        formal_species = {"p": "Proton", "e": "Electron"}

        out = f"{formal_species[self.species]} with energy E = {self.E_eV/1000:.4g}keV."
        return out

    def run(self, orbit=True, info: bool = True, events: list = []):
        r"""
        Calculates the motion and attributes of the particle.

        This functions runs all the required methods for calculating orbit
        in the correct order, and lastly it runs
        :py:meth:`~Particle._orbit`.

        :meta public:

        Parameters
        ----------

        orbit : bool
            Whether or not to actually calculate the orbit in the end. Useful
            when studying particle properties that only depend on initial conditions.
            Defaults to True.
        info : bool
            Whether or not to print the particle's calculated attributes. Defaults to True.
        events : list
            The list of :py:mod:`~gcmotion.scripts.events` to be passed to the solver. Defaults to [].

        """
        logger.info("--------Particle's 'run' routine is called.---------")

        self._conversion_factors()
        self._energies()
        self._orbit_type()

        if orbit:
            start = time()
            solution = self._orbit(events)
            end = time()

            self.theta = solution["theta"]
            self.psi = solution["psi"]
            self.zeta = solution["zeta"]
            self.rho = solution["rho"]
            self.psip = solution["psip"]
            self.Ptheta = solution["Ptheta"]
            self.Pzeta = solution["Pzeta"]
            self.t_eval = solution["t_eval"]
            self.t_events = solution["t_events"]
            self.y_events = solution["y_events"]
            self.message = solution["message"]

            duration = f"{end-start:.4f}"
            self.calculated_orbit = True
            self.solver_output = (
                f"Solver output: {self.message}\n"
                + f"Orbit calculation time: {duration}s."
            )
            logger.info(f"Orbit calculation completed. Took {duration}s")
        else:
            self.solver_output = ""
            logger.info("\tOrbit calculation deliberately skipped.")

        if info:
            logger.info("Printing Particle.__str__() to stdout.")
            print(self.__str__())
        logger.info("Printing Particle.__str__():\n\t\t\t" + self.__str__())

    def _conversion_factors(self):
        r"""
        Calculates the conversion coeffecient needed to convert from lab to NU
        and vice versa.

        Specifically:

        * :math:`\omega_0[s^{-1}] = \dfrac{|Z|e[C]B_0[T]}{m[kg]}`

        * :math:`E_{unit} [J] = m[kg]\omega_0^2[s^{-2}]R^2[m^2]`

        The attributes [NU_to_eV, NU_to_J, Volts_to_NU] are also calculated:

        * To convert from [NU] to [eV], we multiply by
            NU_to_eV = :math:`\dfrac{1}{E_{unit}}`.`

        * To convert from [NU] to [J], we multiply by
            NU_to_J = :math:`\dfrac{e}{E_{unit}}`.`

        * To convert Volts to [NU], we multiply by
            Volts_to_NU = :math:`Z E_{unit}`.`
        """

        logger.info("Calculating conversion factors...")

        qp = physical_constants["qp"]  # 1.6*10**(-19)C
        mp = physical_constants["mp"]  # kg
        B = self.Bfield.B0  # Tesla
        R = self.R  # meters

        self.w0 = qp * B / mp  # [s^-1]
        self.E_unit = mp * self.w0**2 * R**2  # [J]

        # Conversion Factors
        self.NU_to_eV = self.E_unit / qp
        self.NU_to_J = self.E_unit
        self.Volts_to_NU = qp / self.E_unit

        self.calculated_conversion_factors = True
        logger.info("--> Calculated conversion factors.")

    def _energies(self):
        r"""
        Calculates the particle's energy in [NU], [eV] and [J], using
        its initial conditions and the conversion factors.
        """
        r0 = sqrt(2 * self.psi0)
        B_init = self.Bfield.B(r0, self.theta0)
        Phi_init = float(self.Efield.Phi_of_psi(self.psi0))
        Phi_init_NU = Phi_init * self.Volts_to_NU

        self.E_NU = (
            self.qNU**2 / (2 * self.mNU) * self.rho0**2 * B_init**2
            + self.mu * B_init**2
            + self.qNU * Phi_init_NU
        )

        self.E_eV = self.E_NU * self.NU_to_eV
        self.E_keV = self.E_eV / 1000
        self.E_J = self.E_NU * self.NU_to_J

        self.calculated_energies = True
        logger.info("Calculated particle's energies(NU, eV, J).")

    def _orbit_type(self):
        r"""
        Estimates the orbit type given the initial conditions ONLY.

        .. note:: This method works only in the absence of an Electric Field
            and LAR Magnetic Field.

        Trapped/passing:
        The particle is trapped if rho vanishes, so we can
        check if rho changes sign. Since
        :math:`\rho = \dfrac{\sqrt{2W-2\mu B}}{B}`, we need only to
        check under the root.

        Confined/lost:
        (from shape page 87)
        We only have to check if the particle is in-between the 2 left parabolas.
        """
        logger.info("Calculating particle's orbit type:")

        if (self.has_efield) or (not self.Bfield.is_lar):
            self.orbit_type_str = "Cannot calculate (Electric field is present, or Magnetic field is not LAR.)"
            logger.warning(
                "\tElectric field is present, or Magnetic field is not LAR. Orbit type calculation is skipped."
            )
            return

        # Calculate Bmin and Bmax. In LAR, B decreases outwards.
        Bmin = self.Bfield.B(self.a, 0)  # "Bmin occurs at psi_wall, θ = 0"
        Bmax = self.Bfield.B(self.a, np.pi)  # "Bmax occurs at psi_wall, θ = π"

        # Find if trapped or passing from rho (White page 83)
        sqrt1 = 2 * self.E_NU - 2 * self.mu * Bmin
        sqrt2 = 2 * self.E_NU - 2 * self.mu * Bmax
        if sqrt1 * sqrt2 < 0:
            self.t_or_p = "Trapped"
        else:
            self.t_or_p = "Passing"
        logger.debug(f"\tParticle found to be {self.t_or_p}.")

        # Find if lost or confined
        self.orbit_x = self.Pzeta0 / self.psip0
        self.orbit_y = self.mu / self.E_NU
        logger.debug("\tCallling Construct class...")
        foo = Construct(self, get_abcs=True)

        # Recalculate y by reconstructing the parabola (there might be a better way
        # to do this)
        upper_y = (
            foo.abcs[0][0] * self.orbit_x**2
            + foo.abcs[0][1] * self.orbit_x
            + foo.abcs[0][2]
        )
        lower_y = (
            foo.abcs[1][0] * self.orbit_x**2
            + foo.abcs[1][1] * self.orbit_x
            + foo.abcs[1][2]
        )

        if self.orbit_y < upper_y and self.orbit_y > lower_y:
            self.l_or_c = "Confined"
        else:
            self.l_or_c = "Lost"
        logger.debug(f"\tParticle found to be {self.l_or_c}.")

        self.orbit_type_str = self.t_or_p + "-" + self.l_or_c

        self.calculated_orbit_type = True
        logger.info(
            f"--> Orbit type completed. Result: {self.orbit_type_str}."
        )

    def _orbit(self, events: list = []):
        """Groups the particle's initial conditions and passes the to the solver Script
        :py:mod:`~gcmotion.scripts.orbit`. The unpacking takes place in
        :py:mod:`~gcmotion.classes.particle.Particle.run`.

        Parameters
        ----------

        events : list, optional
            List containing the :py:mod:`~gcmotion.scripts.events`. Defaults to []

        Returns
        -------

        dict
            The ``solution`` dictionary returned by the solver
        """

        t = self.t_eval

        init_cond = {
            "theta0": self.theta0,
            "psi0": self.psi0,
            "zeta0": self.zeta0,
            "rho0": self.rho0,
        }

        constants = {
            "mu": self.mu,
            "mi": self.mNU,
            "qi": self.qNU,
        }

        profile = {
            "qfactor": self.qfactor,
            "Bfield": self.Bfield,
            "Efield": self.Efield,
            "Volts_to_NU": self.Volts_to_NU,
        }

        return orbit(t, init_cond, constants, profile, events)
