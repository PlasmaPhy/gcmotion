import numpy as np
import pint
from collections import namedtuple
from time import time

from gcmotion.utils.logger_setup import logger

from gcmotion.utils.pprint_dict import pprint_dict
from gcmotion.utils.get_size import get_size

from gcmotion.scripts.orbit import orbit

from gcmotion.entities.profile import Profile
from gcmotion.entities.initial_conditions import (
    InitialConditions,
    _InitialConditionsFull,
)

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


class Particle(Profile, InitialConditions):
    r"""
    A particle entity represents a fully-fledged particle inside a specific
    tokamak device, and defined initial conditions.

    Parameters
    ----------
    tokamak : Tokamak
        Base object containing information about the tokamak.
    params : PhysicalParameters
        Base object containing the species, :math:`\mu` and :math:`P_\zeta`
        Quantities.
    init : InitialConditions
        Base object containing the :math:`\theta_0, \zeta_0, \psi_0` and
        :math:`t_{eval}` Quantities.

    """

    def __init__(
        self,
        profile: Profile,
        init: InitialConditions,
    ):
        r"""Initializes particle Quantities and Tokamak configuration."""
        logger.info(f"==> Initializing {profile.species_name}...")

        # Grabb attributes from parents
        Profile.__init__(
            self,
            tokamak=profile.tokamak,
            params=profile.params,
        )

        _InitialConditionsFull.__init__(
            self,
            init=init,
            profile=profile,
        )

        # Store those for easier reference
        self.profile = profile
        self.init = init

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
        Calls :py:meth:`~Particle._orbit` to calculate the particle's orbit,
        which returns the solution in [NU]. Then it stores the results in both
        [SI] and [NU].

        :meta public:

        Parameters
        ----------
        orbit : bool, optional
            Whether or not to actually calculate the orbit. Useful when
            studying particle properties that only depend on initial
            conditions. Defaults to True.
        info : bool, optional
            Whether or not to print an output message. Defaults to True.
        events : list, optional
            The list of :py:mod:`~gcmotion.scripts.events` to be passed to the
            solver. Defaults to [].
        """
        logger.info("--------Particle's 'run' routine is called.---------")

        # self._orbit_type()

        if orbit:
            logger.info("Calculating orbit in NU...")
            tfinal = self.t_eval[-1]
            # Orbit Calculation
            start = time()
            solution = self._orbit(events=events)
            end = time()
            solve_time = self.Q(end - start, "seconds")
            logger.info(
                f"Calculation complete. Took {
                        solve_time:.4g~#P}."
            )

            self.theta = self.Q(solution.theta, "radians")
            self.zeta = self.Q(solution.zeta, "radians")
            self.psiNU = self.Q(solution.psi, "NUMagnetic_flux")
            self.rhoNU = self.Q(solution.rho, "NUmeters")
            self.psipNU = self.Q(solution.psip, "NUMagnetic_flux")
            self.PthetaNU = self.Q(solution.Ptheta, "NUMagnetic_flux")
            self.PzetaNU = self.Q(solution.Pzeta, "NUMagnetic_flux")
            self.t_evalNU = self.Q(solution.t_eval, "NUseconds")
            self.t_eventsNU = self.Q(solution.t_events, "NUseconds")
            self.y_events = solution.y_events
            message = solution.message

            # Converting back to SI
            logger.info("Converting results to SI...")
            start = time()
            self.psi = self.psiNU.to("Magnetic_flux")
            self.rho = self.rhoNU.to("meters")
            self.psip = self.psipNU.to("Magnetic_flux")
            self.Ptheta = self.PthetaNU.to("Magnetic_flux")
            self.Pzeta = self.PzetaNU.to("Magnetic_flux")
            self.t_eval = self.t_evalNU.to("seconds")
            self.t_events = self.t_eventsNU.to("seconds")
            end = time()
            conversion_time = self.Q(end - start, "seconds")

            logger.info(
                f"Conversion completed. Took {
                        conversion_time:.4g~#P}."
            )

            # Percentage of t_eval
            percentage = 100 * (self.t_eval[-1] / tfinal).magnitude

            self.solver_output = (
                "\nSolver output:\n"
                + f"{'Message':>23} : "
                + f"{message}\n"
                + f"{'Percentage':>23} : "
                + f"{percentage:.1f}%\n"
                + f"{'Orbit calculation time':>23} : "
                + f"{solve_time:.4g~#P}\n"
                + f"{'Conversion to SI time':>23} : "
                + f"{conversion_time:.4g~#P}\n"
            )

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
        We only have to check if the particle is in-between the 2 left
        parabolas.
        """
        logger.info("Calculating particle's orbit type:")

        if (self.has_efield) or (not self.bfield.is_lar):
            self.orbit_type_str = (
                "Cannot calculate (Electric field is present"
                + "or Magnetic field is not LAR.)"
            )
            logger.warning(
                "\tElectric field is present, or Magnetic field is not LAR. "
                + "Orbit type calculation is skipped."
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

        # # Recalculate y by reconstructing the parabola (there might be a
        # # better way to do this).
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

    def _orbit(self, events: list = []):
        """Groups the particle's initial conditions and passes them to the
        solver Script :py:mod:`~gcmotion.scripts.orbit`. The unpacking takes
        place in :py:meth:`~gcmotion.entities.particle.Particle.run`.

        :meta private:

        Parameters
        ----------
        events : list, optional
            List containing the :py:mod:`~gcmotion.scripts.events`. Defaults to
            []

        Returns
        -------
        namedtuple
            The ``solution`` tuple returned by the solver.
        """
        Parameters = namedtuple(
            "Orbit_Parameters", ["theta0", "psi0", "zeta0", "rho0", "mu", "t"]
        )

        parameters = Parameters(
            theta0=self.theta0,
            psi0=self.psi0NU.magnitude,
            zeta0=self.zeta0,
            rho0=self.rho0NU.magnitude,
            mu=self.muNU.magnitude,
            t=self.t_evalNU.magnitude,
        )

        profile = self.profile

        return orbit(parameters, profile, events=events)

    def __str__(self):
        return "foo"

    def __repr__(self):
        return "foo"

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
