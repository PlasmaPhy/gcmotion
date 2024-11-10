import numpy as np
import pint
from time import time
from tqdm import tqdm

from gcmotion.classes.particle import Particle

# from gcmotion.scripts.events import when_theta_trapped, when_theta_passing

from gcmotion.utils.logger_setup import logger
from gcmotion.utils.get_size import get_size
from gcmotion.scripts.events import when_theta


class Collection:
    r"""
    Instantiates a collection of particles.

    All the particles are stored inside a ``Collection`` object as
    its attributes. We can view them at any time. For example, for a
    collection of 3 Protons:

    .. code-block:: python

        >>> collection.particles #FIXME:
        Proton with energy E=2.723keV.,
        Proton with energy E=2.224keV.,
        Proton with energy E=1.833keV.]

        >>> collection.particles[2]
        Proton with energy E=1.833keV.

        >>> print(collection.particles[2])
        <particle's ``__str__()`` method>

        >>> collection.params
        <parameters imported from ``parameters.py`` file>



    Most of the methods as well as plotters are
    simple wrappers around the
    :py:class:`~gcmotion.classes.particle.Particle` class and the
    single-particle plotters.

    Example
    -------

    Here is a way to initialize a collection of particles and run them:

    .. code-block:: python

        >>> collection = gcm.Collection(parameters)
        >>> collection.run_all(orbit=True, terminal = 10)

    .. rubric:: Methods
        :heading-level: 1
    """

    def __init__(self, file):
        r"""
        Initializes a Collection of particles.

        This class is a *Collection* of many independent particles, each one with
        its own properties, configuration and initial condition. It reads all their
        attributesform ``parameters.py`` file, which must be imported such:

        .. code-block:: python

            >>> from <path-to-file> import parameters

        The file must contain a dict called "params" which contains all the information.
        See example.

        The Collection can then calculate the particle's orbits.

        Parameters
        ----------
        file : .py file
            The python file containing the parameters.

        """

        logger.info("----------------Initializing Collection----------------")

        # Import data and check them
        self.ureg, self.Q = file.ureg, file.Q
        self.params = file.params
        self._check(file)
        self._create()

    def _check(self, file) -> bool:
        r"""
        Checks for the validity of the parameters file.

        Checks if all given parameter lists are of length 1 or otherwise of
        equal length with each other, and if their values are valid.

        Also creates truth flags for every parameter: True if single-valued,
        and False if variable.

        Parameters
        ----------
        file : .py file
            The python file containing the parameters.

        Returns
        -------
        bool
            Truth value depending on the validity of the file

        """

        print("Checking Data...")
        logger.info("Checking Data...")

        # Store the lenghts of each parameter and find the number of particles.
        # The "lengths" dictionary has the same keys as "params", and the values
        # are the length of the corresponding value in "params".
        self.lengths = {x: 1 for x in self.params}
        for key, value in self.params.items():
            if isinstance(value, pint.Quantity) and key != "t_eval":
                if isinstance(value.magnitude, (list, np.ndarray)):
                    self.lengths[key] = len(value.magnitude)
            elif isinstance(value, (list, np.ndarray)):
                self.lengths[key] = len(value)

        # Number of particles "n" is assigned to be the length of the Quantity with
        # the most values. Whether or not every Quantity has the same number of values
        # is checked later.
        self.n = max(self.lengths.values())
        logger.debug(f"Number of particles found: {self.n}")

        # Create a bool flag for every Quantity. True if multiple values are found, and
        # False if single-valued
        for key, value in self.lengths.items():
            setattr(self, "multiple_" + key, bool(value - 1))
            logger.debug(f"Flag: 'multiple_{key}` = {bool(getattr(self, "multiple_" + key))}")  # fmt: skip

        # All Quantities must be either single-valued, or have a total of n values.
        if all((_ == self.n) or (_ == 1) for _ in self.lengths.values()):
            print(f"Data is OK. Number of particles = {self.n}")
            logger.info(f"Data is OK. Number of particles = {self.n}")
            return True
        else:
            print(
                "Error: Multiple valued parameters must all be of same length"
            )
            logger.error("Error: Multiple valued parameters must all be of same length")  # fmt: skip
            return False

    def _create(self):
        r"""Initiates the particles."""

        logger.info("Initializing particles...")

        # Make a local, iterable copy of params by expanding all the single-valued
        # Quantities to repeat n times.
        # CAUTION! all objects of same value point to one single object.
        # Changing one of them changes all of them.
        params = self.params.copy()
        for key, value in params.items():
            if self.lengths[key] == 1:
                params[key] = [value] * self.n

        self.particles = []

        logger.disable("gcmotion")
        for i in range(self.n):
            tokamak = {
                "R": params["R"][i],
                "a": params["a"][i],
                "B0": params["B0"][i],
                "qfactor": params["qfactor"][i],
                "bfield": params["bfield"][i],
                "efield": params["efield"][i],
            }
            parameters = {
                "species": params["species"][i],
                "mu/muB": params["mu/muB"][i],
                "theta0": params["theta0"][i],
                "zeta0": params["zeta0"][i],
                "psi0": params["psi0"][i],
                "Pzeta0": params["Pzeta0"][i],
                "t_eval": params["t_eval"][i],
            }

            # Particle Creation
            p = Particle(tokamak, parameters)

            self.particles.append(p)
        logger.enable("gcmotion")
        logger.info("Particle Initialization complete.")

        # Check if at least one of the particles has an electric field.
        self.has_efield = any((p.has_efield) for p in self.particles)

    def run_all(
        self, orbit=True, terminal: int = 0, pole: int | float = np.pi / 2
    ):
        r"""Calculates all the particle's orbit, by running Particle.run() itself.

        Some plots and calculations, such as the parabolas and the orbit type
        calculation don't require the whole orbit to be calculated, since they
        only depend on the initial conditions. We can therefore save valuable time.

        Also keeps statistics of calculation times and event triggers.

        .. caution::
            The ``terminal`` parameter does not specify the number of *full periods*
            before halting the solver, but rather the number that the :math:`psi`
            coordinate has encountered the same value. For more info, see ``note``
            in :py:func:`~gcmotion.scripts.events.when_psi`.

        Parameters
        ----------

        orbit : bool, optional
            Whether or not to calculate the particles' orbits. Defaults to True.
        terminal : int
            The number of event triggers before stopping the orbit calculation of
            each particle. Defaults to 0, which makes the event non-terminal.
        pole : int | float
            The modulo pole used for the event that stops passing particles.
            Must be different than **ANY** of the initial theta0s, and should
            lie between (0,2π). Defaults to π/2.
        """
        logger.info("Calculating particle's orbits...")

        # Statistics
        times = []
        reached_end = 0
        terminated = 0

        pbar = tqdm(total=self.n)
        for p in self.particles:

            logger.disable("gcmotion")

            start = time()
            events = [
                when_theta(p.theta0, 10)
            ]  # FIXME: update events to work with Quantities
            p.run(info=False, orbit=orbit, events=events)
            times.append(time() - start)

            logger.enable("gcmotion")

            if "0" in p.solver_output:  # Reached the end of integration
                reached_end += 1
            else:  # terminated by the event
                terminated += 1

            pbar.update(1)

        pbar.close()
        logger.enable("gcmotion")

        times = np.array(times)
        time_str = (
            f"\nTotal calculation time: {times.sum():.4g}s.\n"
            + f"Fastest particle: {times.min():.4g}s.\n"
            + f"Slowest particle: {times.max():.4g}s.\n"
            + f"Average time: ({times.mean():.4g} \u00B1 {times.std():.4g})s.\n"
            + f"{reached_end} particles reached the end of the integration interval.\n"
            + f"{terminated} particles were terminated early.\n"
        )
        print(time_str)
        logger.info(time_str)

    def __repr__(self):
        return f"Collection of {self.n} particles."

    def __getitem__(self, item):
        try:
            return self.particles[item]
        except TypeError:
            raise TypeError("Index must be an integer.")

    def __sizeof__(self):
        sizes = [
            self.Q(get_size(vars(self[_])), "bytes") for _ in range(self.n)
        ]
        size = sum(sizes)
        return f"{size:.4g~#P}"
