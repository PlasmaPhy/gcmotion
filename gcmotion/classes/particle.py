import numpy as np
import pint
from time import time
from collections import namedtuple

from gcmotion.tokamak.efield import Nofield

from gcmotion.utils.logger_setup import logger

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
    A particle holds its properties such as mass, species, etc as well as all its calculated
    quantities as its attributes.

    Physical properties and tokamak configuration are automatically setup upon the particle's
    initialization, and its orbit is calculated upon calling its
    :py:meth:`~Particle.run` method.

    Here is a list of the most important attributes:

    #. Initial conditions:
        theta0, psi0, zeta0, rho0, psip0, Ptheta0, Pzeta0
    #. Constants of motion:
        E, EkeV, E, mu, mu
    #. Time evolution arrays:
        theta, psi, zeta, rho, psip, Ptheta, Pzeta
    #. Configuration objects and parameters:
        R, a, B0, qfactor, bfield, efield, psi_wall

    .. hint::
        For every quantitative attribute its corresponing value in [NU] is also
        an attribute which is has the same name as the [SI] attribute with the
        suffix "NU". for example:

        .. code-block:: python

            >>> cwp["Pzeta0"]
            -0.074052 Tm^2
            >>> cwp["Pzeta0NU"]
            -0.0272 NUmf    # NU Magnetic flux


    .. hint::
        We can always view all the currently stored *Quantities* by running:

        .. code-block:: python

            >>> cwp.quantities(everything=True)

    Example
    -------
    See :ref:`Creating a particle <user_guide_particle_creation>`

    .. rubric:: Methods
        :heading-level: 3
    """

    def __init__(
        self,
        tokamak: dict,
        parameters: dict,
    ):
        r"""Initializes particle Quantities and Tokamak configuration.

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
                    [Magnetic flux density].
                qfactor : :py:class:`~gcmotion.tokamak.qfactor.QFactor`
                    Qfactor object that supports query methods for getting values
                    of :math:`q(\psi)` and :math:`\psi_p(\psi)`.
                bfield : :py:class:`~gcmotion.tokamak.bfield.MagneticField`
                    Magnetic Field Object that supports query methods for getting values of the
                    field magnitude, plasma currents and their derivatives.
                efield : :py:class:`~gcmotion.tokamak.efield.ElectricField`
                    Electric Field Object that supports query methods for getting values of the
                    field itself and the derivatives of its potential.
        parameters : dict
            A dict containing all the particle-specific parameters.

                "species" : str
                    The particle species, used to set up charge and mass automatically
                    (from :py:mod:`~gcmotion.configuration.particle_attributes`)
                "mu/muB": Quantity
                    This input parameter is parsed differently depending on its
                    **dimensionality**. If its dimensionality is
                    :math:`[current] [length]^2` (dimensionality of magnetic moment), then
                    it is parsed as the particle's :math:`\mu`. If its dimensionality is
                    :math:`[mass] \dfrac{[length]^2}{[time]^2}` (dimensionality of energy),
                    then it is parsed as the particle's initial :math:`\mu B`. Keep in mind that
                    dimentionality is different than the units system. In both occations, the
                    input can be either in [SI] or [NU], since the dimensionality is
                    independent of the units system used.
                "theta0" and "zeta0" : float, Quantity
                    The :math:`\theta_0, \zeta_0` initial conditions [radians/dimensionless].
                "psi0" : Quantity
                    The :math:`\psi_0` initial values in dimensions of [Magnetic_flux] or
                    [psi_wall]. [psi_wall] is a unit of Magnetic flux, where 1[psi_wall] is
                    defined to be the Magnetic flux of the last closed surface. This way we
                    can set the :math:`\psi_0` initial value with respect to the tokamak's
                    :math:`\psi_{wall}`.
                "Pzeta0" : Quantity
                    The :math:`P_{\zeta_0}` initial value in dimensions of [Magnetic flux].
                "t_eval" : Quantity (np.ndarray)
                    The time interval return values, [:math:`t_0, t_f, steps`], in dimensions
                    of [time].

        """

        logger.info(f"--------Initializing {particle_attributes[parameters["species"]+"_name"]}--------")  # fmt: skip

        def setup_species():
            """Grabs particle's constants from ``particle_attributes.py``
            and sets up its mass and charge."""
            logger.info("Setting up particle's constants...")

            self.species = parameters["species"].lower()
            # Atomic mass and weight, in integer units of proton masses
            M = particle_attributes[self.species + "_M"]  # Purely numeric
            Z = particle_attributes[self.species + "_Z"]  # Purely numeric

            self.miNU = self.Q(M, "Proton_mass")
            self.qiNU = self.Q(Z, "Proton_charge")

            self.mi = self.miNU.to("kilogram")
            self.qi = self.qiNU.to("Coulomb")

        def setup_tokamak():
            """Sets up tokamak configuration and attributes."""
            logger.info("Setting up Tokamak...")

            # Convert and store input in SI, regardless of its initial units
            self.R = tokamak["R"].to("meters")
            self.a = tokamak["a"].to("meters")
            self.B0 = tokamak["B0"].to("Tesla")

            # Store their corresponding NU Quantities
            self.RNU = self.R.to("NUmeters")
            self.aNU = self.a.to("NUmeters")
            self.B0NU = self.B0.to("NUTesla")

            # Objects setup
            self.qfactor = tokamak["qfactor"]
            self.bfield = tokamak["bfield"]
            efield = tokamak["efield"]
            if efield is None or isinstance(efield, Nofield):
                self.efield = Nofield()
                self.has_efield = False
            else:
                self.efield = efield
                self.has_efield = True

            # Last closed surfaces
            # Note: psipNU() input and output are in [NU].
            self.psi_wall = (self.B0 * self.a**2 / 2).to("Magnetic_flux")
            self.psi_wallNU = self.psi_wall.to("NUMagnetic_flux")
            self.psip_wallNU = self.Q(
                self.qfactor.psipNU(self.psi_wallNU.magnitude),
                "NUMagnetic_flux",
            )
            self.psip_wall = self.psip_wallNU.to("Magnetic_flux")

            # Group configuration objects in a named tuple
            # for easier reference
            Profile = namedtuple(
                "Tokamak_Profile", ["qfactor", "bfield", "efield"]
            )
            self.profile = Profile(
                qfactor=self.qfactor,
                bfield=self.bfield,
                efield=self.efield,
            )

        def setup_parameters():
            """Sets up the particles initial conditions, parameters
            and various derivative attributes."""
            logger.info("Setting up particle's initial conditions...")

            # Convert and store input in SI, regardless of its initial units
            self.theta0 = self.Q(parameters["theta0"], "radians")
            self.zeta0 = self.Q(parameters["zeta0"], "radians")
            self.psi0 = self.Q(parameters["psi0"], "Magnetic_flux")
            self.Pzeta0 = parameters["Pzeta0"].to("Magnetic_flux")
            self.t_eval = parameters["t_eval"].to("seconds")
            self.tfinal = self.t_eval[-1]

            # Store their corresponding NU Quantities
            self.psi0NU = self.psi0.to("NUMagnetic_flux")
            self.Pzeta0NU = self.Pzeta0.to("NUmagnetic_flux")
            self.t_evalNU = self.t_eval.to("NUseconds")

            # Initial (B,i,g) and Phi values
            # Note: bigNU() and PhiNU() input and output are in [NU]
            # Symbolize the magnitude of the repspective variable with a _
            _psi0NU = self.psi0NU.magnitude
            _theta0 = self.theta0.magnitude
            B_init = (self.B0 * self.bfield.bigNU(_psi0NU, _theta0)[0]).to(
                "Tesla"
            )
            i_init = self.Q(
                self.bfield.bigNU(_psi0NU, _theta0)[1],
                "NUPlasma_current",
            ).to("Plasma_current")
            g_init = self.Q(
                self.bfield.bigNU(_psi0NU, _theta0)[2],
                "NUPlasma_current",
            ).to("Plasma_current")
            Phi_init = self.Q(
                self.efield.PhiNU(_psi0NU, _theta0),
                "Volts",
            )

            # mu-muB parsing
            # define "mu" and "muB" accordingly depending on which was given by
            # checking the dimensionality of "mu/muB"
            magnetic_moment_dim = self.ureg.Magnetic_moment.dimensionality
            energy_dim = self.ureg.Joule.dimensionality

            if parameters["mu/muB"].dimensionality == magnetic_moment_dim:
                self.mu = parameters["mu/muB"].to("Magnetic_moment")
                self.muB = (self.mu * B_init).to("keV")
            elif parameters["mu/muB"].dimensionality == energy_dim:
                self.muB = parameters["mu/muB"].to("keV")
                self.mu = self.muB / B_init

            self.psip0 = self.Q(  # psipNU only accepts [NU]
                self.qfactor.psipNU(_psi0NU),
                "NUMagnetic_flux",
            ).to("Magnetic_flux")

            self.rho0 = ((self.Pzeta0 + self.psip0) / g_init).to(
                "meters"  # Note: [rho] = Magnetic_flux / Plasma_current = meters
            )

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

        # Grab Q and UnitRegistry
        self.ureg, self.Q = setup_pint(
            R=tokamak["R"].magnitude,
            a=tokamak["a"].magnitude,
            B0=tokamak["bfield"].B0.magnitude,
            species=parameters["species"],
        )

        setup_species()
        setup_tokamak()
        setup_parameters()
        self.input_vars = vars(self).copy()  # Store __init__() vars
        logger.info("--------Particle Initialization Completed--------")

    def quantities(
        self,
        which: str = "",
        everything: bool = False,
    ):
        """Prints the pint Quantities of the object.

        :meta public:

        Parameters
        ----------
        which : str, optional
            Options on which Quantities to print. Can include many options.

                #. "init" :
                    Prints all the Quantities defined upon the particle's
                    instanciation.
                #. "NU" or "SI":
                    Print the respective subset of Quantities

            Options can be chained together, for example "initNU"
            Defaults to "" (prints all *Quantites*)

        everything : bool, optional
            Whether or not to print *all* particle's *attributes*, and not just
            Quantities. Ignores "which" arguement if True. Defaults to False.

        """

        units = "NU" if "NU" in which else "SI" if "SI" in which else ""

        if "init" in which:
            pprint_dict(self.input_vars, everything=everything, units=units)
        else:
            pprint_dict(self.__dict__, everything=everything, units=units)

    def run(
        self,
        orbit=True,
        info: bool = True,
        events: list = [],
    ):
        r"""
        Calls :py:meth:`~Particle._orbit` to calculate the particle's orbit, which
        returns the solution in [NU]. Then it stores the results in both [SI]
        and [NU].

        :meta public:

        Parameters
        ----------
        orbit : bool, optional
            Whether or not to actually calculate the orbit. Useful when studying
            particle properties that only depend on initial conditions.
            Defaults to True.
        info : bool, optional
            Whether or not to print an output message. Defaults to True.
        events : list, optional
            The list of :py:mod:`~gcmotion.scripts.events` to be passed to the solver.
            Defaults to [].
        """
        logger.info("--------Particle's 'run' routine is called.---------")

        # self._orbit_type()

        if orbit:
            logger.info("Calculating orbit in NU...")
            # Orbit Calculation
            start = time()
            solution = self._orbit(events=events)
            end = time()
            solve_time = self.Q(end - start, "seconds")
            logger.info(f"Calculation complete. Took {solve_time:.4g~#P}.")  # fmt: skip

            # fmt: off
            self.theta      = self.Q(solution.theta, "radians")
            self.zeta       = self.Q(solution.zeta, "radians")
            self.psiNU      = self.Q(solution.psi, "NUMagnetic_flux")
            self.rhoNU      = self.Q(solution.rho, "NUmeters")
            self.psipNU     = self.Q(solution.psip, "NUMagnetic_flux")
            self.PthetaNU   = self.Q(solution.Ptheta, "NUMagnetic_flux")
            self.PzetaNU    = self.Q(solution.Pzeta, "NUMagnetic_flux")
            self.t_evalNU   = self.Q(solution.t_eval, "NUseconds")
            self.t_eventsNU = self.Q(solution.t_events, "NUseconds")
            self.y_events   = solution.y_events
            message         = solution.message

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

            # Percentage of t_eval
            percentage = 100*(self.t_eval[-1]/self.tfinal).magnitude

            self.solver_output = (
                "\nSolver output:\n"
                + f"{'Message':>23} : " + f"{message}\n" 
                + f"{'Percentage':>23} : " + f"{percentage:.1f}%\n"
                + f"{'Orbit calculation time':>23} : " + f"{solve_time:.4g~#P}\n"
                + f"{'Conversion to SI time':>23} : " + f"{conversion_time:.4g~#P}\n"
            )
            # fmt: on
        else:
            self.solver_output = (
                "\nSolver output: Orbit Calculation deliberately skipped.\n"
            )
            logger.info("\tOrbit calculation deliberately skipped.")

        if info:
            print(self.__str__())

        logger.info("--------Particle's 'run' routine returned.---------")

    def _orbit_type(self):  # FIXME: Parabolas need re-writing.
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

    def _orbit(self, events: list = []) -> namedtuple:
        """Groups the particle's initial conditions and passes them to the solver Script
        :py:mod:`~gcmotion.scripts.orbit`. The unpacking takes place in
        :py:meth:`~gcmotion.classes.particle.Particle.run`.

        :meta private:

        Parameters
        ----------
        events : list, optional
            List containing the :py:mod:`~gcmotion.scripts.events`. Defaults to []

        Returns
        -------
        namedtuple
            The ``solution`` tuple returned by the solver.
        """
        Parameters = namedtuple(
            "Orbit_Parameters", ["theta0", "psi0", "zeta0", "rho0", "mu", "t"]
        )

        # fmt: off
        parameters = Parameters(
            theta0 = self.theta0.magnitude,
            psi0   = self.psi0NU.magnitude,            
            zeta0  = self.zeta0.magnitude,
            rho0   = self.rho0NU.magnitude,
            mu     = self.muNU.magnitude,
            t      = self.t_evalNU.magnitude
        )
        # fmt: on

        profile = self.profile

        return orbit(parameters, profile, events=events)

    def __str__(self):
        delimeter = "\n" + "=" * 100 + "\n"
        particle_name = particle_attributes[self.species + "_name"]
        solver_output = getattr(self, "solver_output", "")

        # fmt: off
        tokamak = (
            "\nTokamak:\n"
            + f"{'R':>23} : " + f"{f'{self.R:.4g~P}':<16}" + f"({self.RNU:.4g~P})"  + "\n"
            + f"{'a':>23} : " + f"{f'{self.a:.4g~P}':<16}" + f"({self.aNU:.4g~P})" + "\n"
            + f"{'B0':>23} : " + f"{f'{self.B0:.4g~P}':<16}" + f"({self.B0NU:.4g~P})" + "\n"
            + f"{'q-factor':>23} : " + f"{self.qfactor}" + "\n"
            + f"{'Magnetic Field':>23} : " + f"{self.bfield}" + "\n"
            + f"{'Electric Field':>23} : " + f"{self.efield}" + "\n"
        )

        particle = (
            "\nParticle:\n"
            + f"{'Particle species':>23} : " + particle_name +"\n"
            + f"{'μB product':>23} : " + f"{f'{self.muB:.4g~P}':<16}" + f"({self.muBNU:.4g~P})" + "\n"
            + f"{'Mangetic moment':>23} : " + f"{f'{self.mu:.4g~P}':<16}" + f"({self.muNU:.4g~P})" + "\n"
            + f"{'Pzeta':>23} : " + f"{f'{self.Pzeta0:.4g~P}':<16}" + f"({self.Pzeta0NU:.4g~P})" + "\n"
            + f"{'Energy (keV)':>23} : " + f"{f'{self.EkeV:.4g~P}':<16}" + f"({self.ENU:.4g~P})" + "\n"
        )
        # fmt: on

        return delimeter + tokamak + particle + solver_output + delimeter

    def __repr__(self):
        particle_name = particle_attributes[self.species + "_name"]
        return f"{particle_name}: E={self.EkeV:.4g~P}, muB={self.muB:.4g~P}, Pzeta={self.Pzeta0:.4g~P}"

    def __getitem__(self, item):
        item = getattr(self, item)
        if isinstance(item, pint.Quantity):
            print(f"{item:.6g~}")
        else:
            return item

    def __sizeof__(self):
        r"""Recursively calculates the size of the instanciated particle.

        Might take a couple of seconds.

        :meta public:

        """

        size = self.Q(get_size(vars(self)), "bytes")
        return f"{size:.4g~#P}"
