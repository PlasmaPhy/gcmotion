import numpy as np
import pint
from time import time
from pprint import pprint

from gcmotion.tokamak.efield import Nofield

from gcmotion.utils._logger_setup import logger

from gcmotion.configuration.particle_attributes import particle_attributes
from gcmotion.utils.setup_pint import setup_pint

from gcmotion.scripts.orbit import orbit

# from gcmotion.classes.parabolas import Construct


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
                (from :py:mod:`~gcmotion.configuration.particle_attributes`)

        """
        logger.info("--------Initializing particle--------")

        def setup_constants():
            """Grabs particle's constants from ``particle_attributes.py``"""

            logger.info("Setting up particle's constants...")

            self.species = parameters["species"].lower()

            M = particle_attributes[self.species + "_mNU"]
            Z = particle_attributes[self.species + "_qNU"]

            self.mNU = self.Q(M, "Proton_mass")
            self.qNU = self.Q(Z, "Proton_charge")

            self.mi = self.mNU.to("kilogram")
            self.qi = self.qNU.to("Coulomb")

            logger.debug(f"\tParticle is of species '{self.species}'.")
            logger.info("--> Particle's constants setup successful")

        def setup_tokamak():
            """Sets up tokamak-related attributes."""

            logger.info("Setting up Tokamak...")

            # Dimensions
            self.R = self.Q(tokamak["R"], "meters")
            self.a = self.Q(tokamak["a"], "meters")

            logger.debug(f"\tTokamak dimensions: R = {self.R}, a = {self.a}")

            # Objects
            self.qfactor = tokamak["qfactor"]

            self.Bfield = tokamak["Bfield"]
            self.B0 = self.Q(self.Bfield.B0, "Tesla")
            self.g = self.Q(self.Bfield.g, "Tesla * m")
            self.i = self.Q(self.Bfield.i, "Tesla * m")
            self.B0NU = self.Q(
                self.mNU.magnitude / self.qNU.magnitude, "NUTesla"
            )

            Efield = tokamak["Efield"]
            if Efield is None or isinstance(Efield, Nofield):
                self.Efield = Nofield()
                self.has_efield = False
            else:
                self.Efield = Efield
                self.has_efield = True

            logger.debug(f"\t'{self.qfactor.id}' qfactor used with parameters {self.qfactor.params}")  # fmt: skip
            logger.debug(f"\t'{self.Bfield.id}' Bfield used with parameters {self.Bfield.params}")  # fmt: skip
            logger.debug(f"\t'{self.Efield.id}' Efield used with parameters {self.Efield.params}")  # fmt: skip

            self.psi_wall = self.B0 * self.a**2 / 2  # Tesla * m^2
            self.psip_wall = self.Q(
                self.qfactor.psip_of_psi(self.psi_wall.magnitude),
                self.psi_wall.units,
            )  # Tesla * m^2

            logger.info("--> Tokamak setup successful.")

        def setup_parameters():
            """Sets up the particles initial condition and parameters, as well as the solver's S0."""

            logger.info("Setting up particle's initial conditions...")

            muB = self.Q(parameters["muB"], "keV")
            self.mu = muB / self.B0  # keV / Tesla

            self.theta0 = self.Q(parameters["theta0"], "radians")
            self.zeta0 = self.Q(parameters["zeta0"], "radians")

            self.Pzeta0 = self.Q(
                parameters["Pzeta0"], "NUCanonical_momentum"
            ).to("Tesla * meter ** 2")

            # CAUTION! psi was difined with respect to psi_wall
            self.psi0 = parameters["psi0"] * self.psi_wall  # Tesla * m^2
            self.psip0 = self.Q(
                self.qfactor.psip_of_psi(self.psi0.magnitude),
                self.psi0.units,
            )  # Tesla * m^2

            # [rho] = Magnetic_flux / Plasma_current = meters
            self.rho0 = (self.Pzeta0 + self.psip0) / self.g
            self.rho0 = self.rho0.to_base_units()

            self.Ptheta0 = self.psi0 + self.rho0 * self.i  # Tesla * meter ** 2

            self.t_eval = self.Q(parameters["t_eval"], "seconds")

            logger.debug("ODE initial conditions:\n"+ f"\ttheta0 = {self.theta0:.5g}, psi0 = {self.psi0:.5g}, zeta0 = {self.zeta0:.5g}, Pzeta0 = {self.Pzeta0:.5g}.")  # fmt: skip
            logger.debug("\tOther initial conditions:\n"+ f"\tPtheta0 = {self.Ptheta0:.5g}, psip0 = {self.psip0:.5g}, rho0 = {self.rho0:.5g}, mu = {self.mu:.5g}")  # fmt: skip
            logger.debug(f"\tTime span (t0, tf, steps): ({self.t_eval[0]}, {self.t_eval[-1]}, {len(self.t_eval)})")  # fmt: skip
            logger.info("--> Initial conditions setup successful.")

        def setup_logic_flags():
            """Sets up logic flags and initializes variables that must have an initial value"""

            logger.info("Setting up logic flags...")

            self.t_or_p = "Unknown"
            self.l_or_c = "Unknown"

            # Stored initially to avoid attribute errors
            self.z_0freq = self.z_freq = self.theta_0freq = self.theta_freq = (
                None
            )

            logger.info("--> Logic flags setup successful.")

        _, self.Q = setup_pint(tokamak["R"], tokamak["Bfield"].B0)
        setup_constants()
        setup_tokamak()
        setup_parameters()
        setup_logic_flags()

        logger.info("--------Particle Initialization Completed--------\n")

    def __str__(self):
        # info_str = (
        #     "Constants of motion:\n"
        #     + "\tParticle Energy (normalized):\tE = {:e}\n".format(self.ENU)
        #     + "\tParticle Energy (eV):\t\tE = {:e} eV\n".format(self.E_eV)
        #     + "\tParticle Energy (J):\t\tE = {:e} J\n".format(self.E_J)
        #     + f"\tToroidal Momenta:\t\tPζ = {self.Pzeta0}\n\n"
        #     + "Other Quantities:\n"
        #     + f'\tParticle of Species:\t\t"{self.species}"\n'
        #     + f"\tOrbit Type:\t\t\t{self.orbit_type_str}\n"
        #     + f"\tMajor Radius:\t\t\tR = {self.R} meters\n"
        #     + f"\tMinor Radius:\t\t\tα = {self.a} meters\n"
        #     + "\tToroidal Flux at wall:\t\tψ = {:n}\n".format(self.psi_wall)
        #     + "\tTime unit:\t\t\tω = {:e} Hz \n".format(self.w0)
        #     + "\tEnergy unit:\t\t\tE = {:e} J \n\n".format(self.E_unit)
        #     + self.solver_output
        # )

        info_str = "foo"

        return info_str

    def __repr__(self):
        formal_species = {"p": "Proton", "e": "Electron"}

        out = f"{formal_species[self.species]} with energy E = {self.E_eV/1000:.4g}keV."
        return out

    def quantities(self):
        """Prints the pint Quantities of the object"""
        d = self.__dict__.copy()
        del d["t_eval"]
        d["t_eval_start"] = self.t_eval[:5]
        d["t_eval_end"] = self.t_eval[-5:]
        pprint(
            [
                f"{key} = {value:.4g~}"
                for key, value in d.items()
                if isinstance(value, pint.Quantity)
            ]
        )

    def run(
        self,
        orbit=True,
        units: str = "SI",
        info: bool = True,
        events: list = [],
    ):
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

        # self._conversion_factors()
        self._energies()
        self._orbit_type()

        if orbit:
            start = time()
            solution = self._orbit(events=events, units=units)
            end = time()

            self.theta = solution["theta"]
            self.psi = solution["psi"]
            self.zeta = solution["zeta"]
            self.rho = solution["rho"]
            self.psip = solution["psip"]
            self.Ptheta = solution["Ptheta"]
            self.Pzeta = solution["Pzeta"]
            self.t_eval_sol = solution["t_eval"]  # Returns them dimensionless
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

    def _energies(self):
        r"""
        Calculates the particle's energy in [NU], [eV] and [J], using
        its initial conditions and the conversion factors.
        """
        r0 = (2 * self.psi0 / self.B0) ** (1 / 2)  # math.sqrt doesn't work

        B_init = self.B0 * self.Bfield.B(r0.magnitude, self.theta0.magnitude)
        Phi_init = self.Q(self.Efield.Phi_of_psi(self.psi0.magnitude), "Volts")
        self.E = (
            self.qi**2 / (2 * self.mi) * self.rho0**2 * B_init**2
            + self.mu * B_init
            + self.qi * Phi_init
        ).to("Joule")

        self.EeV = self.E.to("eV")
        self.EkeV = self.E.to("keV")
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
        # "Bmin occurs at psi_wall, θ = 0"
        Bmin = self.B0 * self.Bfield.B(self.a.magnitude, 0)
        # "Bmax occurs at psi_wall, θ = π"
        Bmax = self.B0 * self.Bfield.B(self.a.magnitude, np.pi)

        # Find if trapped or passing from rho (White page 83)
        sqrt1 = 2 * self.E - 2 * self.mu * Bmin
        sqrt2 = 2 * self.E - 2 * self.mu * Bmax
        if sqrt1 * sqrt2 < 0:
            self.t_or_p = "Trapped"
        else:
            self.t_or_p = "Passing"
        logger.debug(f"\tParticle found to be {self.t_or_p}.")

        # # Find if lost or confined
        # self.orbit_x = self.Pzeta0 / self.psip0
        # self.orbit_y = self.mu / self.E
        # logger.debug("\tCallling Construct class...")
        # foo = Construct(self, get_abcs=True)

        # # Recalculate y by reconstructing the parabola (there might be a better way
        # # to do this)
        # upper_y = (
        #     foo.abcs[0][0] * self.orbit_x**2
        #     + foo.abcs[0][1] * self.orbit_x
        #     + foo.abcs[0][2]
        # )
        # lower_y = (
        #     foo.abcs[1][0] * self.orbit_x**2
        #     + foo.abcs[1][1] * self.orbit_x
        #     + foo.abcs[1][2]
        # )

        # if self.orbit_y < upper_y and self.orbit_y > lower_y:
        #     self.l_or_c = "Confined"
        # else:
        #     self.l_or_c = "Lost"
        # logger.debug(f"\tParticle found to be {self.l_or_c}.")

        self.orbit_type_str = self.t_or_p + "-" + self.l_or_c

        self.calculated_orbit_type = True
        logger.info(
            f"--> Orbit type completed. Result: {self.orbit_type_str}."
        )

    def _orbit(self, events: list = [], units: str = "NU"):
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

        t = self.t_eval.to("NUsecond")

        # fmt: off
        parameters = {
            "theta0" : self.theta0,
            "zeta0"  : self.zeta0,
            "psi0"   : self.psi0.to("NUMagnetic_flux"),            
            "rho0"   : self.rho0.to("NUmeter"),
            "mu"     : self.mu.to("NUMagnetic_moment"),
            "mi"     : self.mi.to("Proton_mass"),
            "qi"     : self.qi.to("Proton_charge"),
            "B0"     : self.B0.to("NUTesla"),
            "VtoVNU" : self.Q("1V").to("NUVolts"),
        }

        profile = {
            "qfactor": self.qfactor,
            "Bfield": self.Bfield,
            "Efield": self.Efield,
        }
        # fmt: on

        pprint(
            [
                f"{key} = {value:.4g~}"
                for key, value in parameters.items()
                if isinstance(value, pint.Quantity)
            ]
        )

        return orbit(t, parameters, profile, units=units, events=events)
