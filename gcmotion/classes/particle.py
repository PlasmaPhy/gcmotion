import numpy as np
import pint
from time import time

from gcmotion.tokamak.efield import Nofield

from gcmotion.utils._logger_setup import logger

from gcmotion.configuration.particle_attributes import particle_attributes
from gcmotion.utils.setup_pint import setup_pint
from gcmotion.utils.pprint_dict import pprint_dict
from gcmotion.utils.get_size import get_size

from gcmotion.scripts.orbit import orbit

# from gcmotion.classes.parabolas import Construct

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


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
        R, a, qfactor, bfield, efield, psi_wall
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

                R : Quantity
                    The tokamak's major radius dimensions of [length].
                a : Quantity
                    The tokamak's minor radius dimensions of [length].
                B0 : Quantity
                    The Magnetic field strength on the magnetic axis in dimensions of
                    [Magnetic field strength].
                qfactor : :py:class:`~gcmotion.tokamak.qfactor.QFactor`
                    Qfactor object that supports query methods for getting values
                    of :math:`q(\psi)` and :math:`\psi_p(\psi)`.
                bfield : :py:class:`~gcmotion.tokamak.bfield.MagneticField`
                    Magnetic Field Object that supports query methods for getting values of the
                    field magnitude and its derivatives.
                efield : :py:class:`~gcmotion.tokamak.efield.ElectricField`
                    Electric Field Object that supports query methods for getting values of the
                    field itself and the derivatives of its potential.
        parameters : dict
            A dict containing all the particle-specific parameters.

                species : str
                    The particle species, used to later set charge and mass automatically
                    (from :py:mod:`~gcmotion.configuration.particle_attributes`)
                mu : Quantity
                    The magnetic moment in dimensions of [magnetic moment].
                theta0, zeta0 : float
                    The :math:`\theta_0, \zeta_0` initial conditions (radians-dimensionless).
                psi0, Pzeta0 : Quantity
                    The :math:`\psi_0, P_{\zeta_0}` initial conditons in dimensions of [Magnetic flux].
                t_eval : Quantity (np.ndarray)
                    The ODE time interval return values, [:math:`t_0, t_f`, steps], in dimensions
                    of [time].

        """
        logger.info(f"--------Initializing {particle_attributes[parameters["species"]+"_name"]}--------")  # fmt: skip

        def setup_input():
            r"""Parses input and sets all parameters  as attributes in both SI and NU units"""
            logger.info("Setting up input parameters in SI and NU...")
            # SI attributes
            self.R = tokamak["R"].to("meters")
            self.a = tokamak["a"].to("meters")
            self.B0 = tokamak["B0"].to("Tesla")
            self.theta0 = parameters["theta0"]
            self.zeta0 = parameters["zeta0"]

            # CAUTION! psi was defined with respect to psi_wall
            self.psi_wall = (self.B0 * self.a**2 / 2).to("Magnetic_flux")
            self.psi0 = parameters["psi0"] * self.psi_wall
            self.Pzeta0 = parameters["Pzeta0"].to("Magnetic_flux")
            self.t_eval = parameters["t_eval"].to("seconds")

            # Corresponding NU attributes
            self.RNU = self.R.to("NUmeters")
            self.aNU = self.a.to("NUmeters")
            self.B0NU = self.B0.to("NUTesla")

            self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")
            self.psi0NU = self.psi0.to("NUMagnetic_flux")
            self.Pzeta0NU = self.Pzeta0.to("NUmagnetic_flux")
            self.t_evalNU = self.t_eval.to("NUseconds")

        def setup_species():
            """Grabs particle's constants from ``particle_attributes.py``"""
            logger.info("Setting up particle's constants...")

            self.species = parameters["species"].lower()
            M = particle_attributes[self.species + "_M"]  # Purely numeric
            Z = particle_attributes[self.species + "_Z"]  # Purely numeric

            self.miNU = self.Q(M, "Proton_mass")
            self.qiNU = self.Q(Z, "Proton_charge")

            self.mi = self.miNU.to("kilogram")
            self.qi = self.qiNU.to("Coulomb")

        def setup_tokamak():
            """Sets up tokamak-related attributes."""
            logger.info("Setting up Tokamak...")

            self.qfactor = tokamak["qfactor"]
            self.bfield = tokamak["bfield"]

            efield = tokamak["efield"]
            if efield is None or isinstance(efield, Nofield):
                self.efield = Nofield()
                self.has_efield = False
            else:
                self.efield = efield
                self.has_efield = True

            self.psip_wall = self.Q(
                self.qfactor.psip_of_psi(self.psi_wall.magnitude),
                self.psi_wall.units,
            )
            self.psip_wallNU = self.psip_wall.to("NUmagnetic_flux")

        def setup_parameters():
            """Sets up the particles initial condition and parameters, as well as the solver's S0."""
            logger.info("Setting up particle's initial conditions...")

            # Mu parsing
            r0 = ((self.psi0 / self.B0) ** (1 / 2)).to("meters")
            B_init = (self.B0 * self.bfield.b(r0.magnitude, self.theta0)).to(
                "Tesla"
            )
            if "mu" in parameters.keys():
                self.mu = parameters["mu"].to("Magnetic_moment")
                self.muB = (self.mu * B_init).to("keV")
            elif "muB" in parameters.keys():
                self.muB = parameters["muB"].to("keV")
                self.mu = self.muB / B_init

            # SI attributes
            self.psip0 = self.Q(
                self.qfactor.psip_of_psi(self.psi0NU.magnitude),
                self.psi0NU.units,
            ).to("Magnetic_flux")

            self.rho0 = ((self.Pzeta0 + self.psip0) / self.bfield.g).to(
                "meters"  # Note: [rho] = Magnetic_flux / Plasma_current = meters
            )

            self.Ptheta0 = (self.psi0 + self.rho0 * self.bfield.i).to(
                "Magnetic_flux"
            )

            # Corresponding NU attributes
            self.muNU = self.mu.to("NUMagnetic_moment")
            self.muBNU = self.muB.to("NUkeV")
            self.psip0NU = self.psip0.to("NUMagnetic_flux")
            self.rho0NU = self.rho0.to("NUMeters")
            self.Ptheta0NU = self.Ptheta0.to("NUMagnetic_flux")

        self.ureg, self.Q = setup_pint(
            R=tokamak["R"].magnitude,
            B0=tokamak["bfield"].B0.magnitude,
            species=parameters["species"],
        )

        setup_input()
        setup_species()
        self.input_vars = vars(self).copy()  # Store input values
        setup_tokamak()
        setup_parameters()
        self.clean_size = self.Q(get_size(self), "bytes")

        logger.info("--------Particle Initialization Completed--------\n")

    def __str__(self):
        particle_name = particle_attributes[self.species + "_name"]
        size = self.Q(get_size(self), "bytes")
        EkeV = f"{self.EkeV:.4g~#P}" if hasattr(self, "EkeV") else "Unknown"

        solver_output = getattr(self, "solver_output", "")

        # fmt: off
        tokamak = (
            "Tokamak:\n"
            + f"{'R':>25} : " + f"{self.R:.4g~P}" + "\n"
            + f"{'a':>25} : " + f"{self.a:.4g~P}" + "\n"
            + f"{'B0':>25} : " + f"{self.B0:.4g~P}" + "\n"
            + f"{'q-factor':>25} : " + f"{self.qfactor}" + "\n"
            + f"{'Magnetic Field':>25} : " + f"{self.bfield}" + "\n"
            + f"{'Electric Field':>25} : " + f"{self.efield}" + "\n"
        )

        particle = (
            "Particle:\n"
            + f"{'Particle species':>25} : " + particle_name +"\n"
            + f"{'Energy (keV)':>25} : "+ f"{EkeV}\n"
        )

        diagnostics = (
            "\nDiagnostics:\n"
            + f"{'Init particle size':>25} : " + f"{self.clean_size:.3g~#P}\n"
            + f"{'Full particle size':>25} : " + f"{size:.3g~#P}\n"
        )
        # fmt: on

        return tokamak + particle + diagnostics + solver_output

    def quantities(
        self,
        which="",
        everything=False,
    ):
        """Prints the pint Quantities of the object"""

        units = "NU" if "NU" in which else "SI" if "SI" in which else ""

        if "input" in which:
            pprint_dict(self.input_vars, everything=everything, units=units)
        else:
            pprint_dict(self.__dict__, everything=everything, units=units)

    def run(
        self,
        orbit=True,
        info: bool = False,
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
            Whether or not to print the particle's calculated attributes. Defaults to False.
        events : list
            The list of :py:mod:`~gcmotion.scripts.events` to be passed to the solver. Defaults to [].

        """
        logger.info("--------Particle's 'run' routine is called.---------")

        self._energies()
        # self._orbit_type()

        if orbit:
            logger.info("Calculating orbit in NU...")
            start = time()
            solution = self._orbit(events=events)
            end = time()
            solve_time = self.Q(end - start, "seconds")
            logger.info(f"Calculation complete. Took {solve_time:.4g~#P}.")  # fmt: skip

            # fmt: off
            self.theta      = self.Q(solution["theta"], "radians")
            self.zeta       = self.Q(solution["zeta"], "radians")
            self.psiNU      = self.Q(solution["psi"], "NUMagnetic_flux")
            self.rhoNU      = self.Q(solution["rho"], "NUmeters")
            self.psipNU     = self.Q(solution["psip"], "NUMagnetic_flux")
            self.PthetaNU   = self.Q(solution["Ptheta"], "NUMagnetic_flux")
            self.PzetaNU    = self.Q(solution["Pzeta"], "NUMagnetic_flux")
            self.t_evalNU   = self.Q(solution["t_eval"], "NUseconds")
            self.t_eventsNU = self.Q(solution["t_events"], "NUseconds")
            self.y_events   = solution["y_events"]
            message         = solution["message"]

            # Converting back to SI
            logger.info("Converting results to SI...")
            start = time()
            self.psi        = self.psiNU.to("Magnetic_flux")
            self.rho        = self.rhoNU.to("meters")
            self.psip       = self.psipNU.to("Magnetic_flux")
            self.Ptheta     = self.PthetaNU.to("Magnetic_flux")
            self.Pzeta      = self.PzetaNU.to("Magnetic_flux")
            self.t_eval     = self.t_evalNU.to("seconds")
            self.t_events   = self.t_eventsNU.to("seconds")
            end = time()
            conversion_time = self.Q(end-start, "seconds")

            logger.info(f"Conversion completed. Took {conversion_time:.4g~#P}.")  # fmt: skip

            self.solver_output = (
                "\nSolver output:\n"
                + f"{'Message':>25} : " + f"{message}\n" 
                + f"{'Orbit calculation time':>25} : " + f"{solve_time:.4g~#P}\n"
                + f"{'Conversion to SI time':>25} : " + f"{conversion_time:.4g~#P}"
            )
            # fmt: on
        else:
            self.solver_output = (
                "\nSolver output: Orbit Calculation deliberately skipped.\n"
            )
            logger.info("\tOrbit calculation deliberately skipped.")

        if info:
            print(self.__str__())

    def _energies(self):
        r"""
        Calculates the particle's energy in [NU], [eV] and [J], using
        its initial conditions and the conversion factors.
        """
        logger.info("Calculating Energies...")
        # SI Energies
        r_init = (2 * self.psi0 / self.B0) ** (1 / 2)

        b_init = self.bfield.b(r_init.magnitude, self.theta0)
        B_init = self.B0 * b_init
        Phi_init = self.Q(self.efield.Phi_of_psi(self.psi0.magnitude), "Volts")

        self.E = (
            self.qi**2 / (2 * self.mi) * self.rho0**2 * B_init**2
            + self.mu * B_init
            + self.qi * Phi_init
        ).to("Joule")

        self.EeV = self.E.to("eV")
        self.EkeV = self.E.to("keV")

        # Corresponding NU Energy
        self.ENU = self.E.to("NUJoule")

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

        if (self.has_efield) or (not self.bfield.is_lar):
            self.orbit_type_str = "Cannot calculate (Electric field is present, or Magnetic field is not LAR.)"
            logger.warning(
                "\tElectric field is present, or Magnetic field is not LAR. Orbit type calculation is skipped."
            )
            return

        # Calculate Bmin and Bmax. In LAR, B decreases outwards.
        # "Bmin occurs at psi_wall, θ = 0"
        Bmin = self.B0 * self.bfield.B(self.a.magnitude, 0)
        # "Bmax occurs at psi_wall, θ = π"
        Bmax = self.B0 * self.bfield.B(self.a.magnitude, np.pi)

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

    def _orbit(self, events: list = [], units: str = "SI"):
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
        }

        profile = {
            "qfactor": self.qfactor,
            "bfield" : self.bfield,
            "efield" : self.efield,
        }
        # fmt: on

        return orbit(t, parameters, profile, events=events)
