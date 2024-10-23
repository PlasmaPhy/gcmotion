import numpy as np
from time import time
from tqdm import tqdm

from gcmotion.classes.particle import Particle
from gcmotion.scripts.events import when_theta

from gcmotion.utils._logger_setup import logger


class Collection:
    r"""
    Instantiates a collection of particles.

    All the particles are stored inside a ``Collection`` object as
    its attributes. We can view them at any time. For example, for a
    collection of 3 Protons:

    .. code-block:: python

        >>> collection.particles
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

        # Check data
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
        self.lengths = {x: 1 for x in self.params}
        for key, value in self.params.items():
            if isinstance(value, (int, float)) or key == "t_eval":
                continue
            if isinstance(value, (list, np.ndarray)):
                self.lengths[key] = len(value)

        self.n = max(self.lengths.values())
        logger.debug(f"Number of particles found: {self.n}")

        # Check for multiple efields, bfields, qs, species and create flags
        for key, value in self.lengths.items():
            exec("self.multiple_" + key + "=bool(" + str(value) + "-1)")
            logger.debug(f"Flag: 'multiple_'{key} = {bool(value-1)}")

        # Check lengths and print results
        if all((_ == self.n) or (_ == 1) for _ in self.lengths.values()):
            print(f"Data is OK. Number of particles = {self.n}")
            logger.info(f"Data is OK. Number of particles = {self.n}")
            return True
        else:
            print(
                "Error: Multiple valued parameters must all be of same length"
            )
            logger.error(
                "Error: Multiple valued parameters must all be of same length"
            )
            return False

    def _create(self):
        r"""Initiates the particles."""

        logger.info("Initializing particles...")

        # Make an iterable copy of params
        # CAUTION! all objects of same value point to one single object.
        # Changing one of them changes all of them.
        params = self.params.copy()
        for key, value in params.items():
            if self.lengths[key] == 1:
                params[key] = [value] * self.n

        self.particles = []

        logger.disable("gcmotion")
        for i in range(self.n):
            R, a = params["R"][i], params["a"][i]  # Major/Minor Radius in [m]
            q = params["q"][i]
            Bfield = params["Bfield"][i]
            Efield = params["Efield"][i]

            # Create Particle
            species = params["species"][i]
            mu = params["mu"][i]  # Magnetic moment
            theta0 = params["theta0"][i]
            psi0 = params["psi0"][i]  # times psi_wall
            zeta0 = params["z0"][i]
            Pzeta0 = params["Pz0"][i]
            t_eval = params["t_eval"][i]  # t0, tf, steps

            tokamak = {
                "R": R,
                "a": a,
                "q": q,
                "Bfield": Bfield,
                "Efield": Efield,
            }
            init_cond = {
                "theta0": theta0,
                "psi0": psi0,
                "zeta0": zeta0,
                "Pzeta0": Pzeta0,
            }

            # Particle Creation
            try:
                p = Particle(tokamak, t_eval, init_cond, mu, species)
            except:  # noqa: E722
                error_str = (
                    f"Error initialzing Particle #{i+1}. Dumping:\n"
                    + f"Tokamak condifuration: {tokamak}\n"
                    + f"t_eval(t0, tf, step): {t_eval[0], t_eval[-1], t_eval[1]-t_eval[0]}\n"
                    + f"Initial conditions: {init_cond}\n"
                    + f"mu, species: {mu, species}\n"
                )
                print(error_str)
                logger.enable("gcmotion")
                logger.error(error_str)

            self.particles.append(p)
        logger.enable("gcmotion")
        logger.info("Particle Initialization complete.")

    def _check_multiples(self, allowed: list) -> bool:
        r"""Checks if the given parameters are static or vary from particle to particle.

        Since all the plots require some parameters to be the same and allowing only
        certain parameters to vary from particle to particle, this check is important.
        Otherwise the resulting plots are nonsense.

        Parameters
        ----------
        allowed : list
            The parameters which are alloed to vary.

        Returns
        --------
        bool
            Truth value depending on the result of the check.
        """
        for key in allowed:
            expr = "self.multiple_" + key + " is True"
            if eval(expr):
                print("Error")
                return False
        return True

    def run_all(self, orbit=True, terminal=0):
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
            event = when_theta(p.theta0, terminal)
            p.run(info=False, orbit=orbit, events=[event])
            times.append(time() - start)

            logger.enable("gcmotion")

            if p.message[0] == "0":  # Reached the end of termination integral
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
