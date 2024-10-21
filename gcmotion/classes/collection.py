import numpy as np
from time import time
from tqdm import tqdm

from gcmotion.classes.particle import Particle
from gcmotion.scripts.events import when_psi

from gcmotion.utils._logger_setup import logger, add_logger


class Collection:

    def __init__(self, file):
        """
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
        """
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
            print("Error: Multiple valued parameters must all be of same length")
            logger.error("Error: Multiple valued parameters must all be of same length")
            return False

    def _create(self):
        """Initiates the particles."""

        logger.info("Initializing particles...")

        # Make an iterable copy of params
        # CAUTION! all objects of same value point to one single object.
        # Changing one of them changes all of them.
        params = self.params.copy()
        for key, value in params.items():
            if self.lengths[key] == 1:
                params[key] = [value] * self.n

        self.particles = []

        logger.remove()
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

            tokamak = {"R": R, "a": a, "q": q, "Bfield": Bfield, "Efield": Efield}
            init_cond = {"theta0": theta0, "psi0": psi0, "zeta0": zeta0, "Pzeta0": Pzeta0}

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
                add_logger()
                logger.error(error_str)

            self.particles.append(p)
        add_logger()
        logger.info("Particle Initialization complete.")

    def run_all(self, orbit=True, terminal=0):
        """Calculates all the particle's orbit, by running Particle.run() itself.

        Some plots and calculations, such as the parabolas and the orbit type
        calculation don't require the whole orbit to be calculated, since they
        only depend on the initial conditions. We can therefore save valuable time.

        Also keeps statistics of calculation times and event triggers.

        Parameters
        ----------

        orbit : bool, optional
            Whether or not to calculate the particles' orbits. Defaults to True.
        terminal : int
            The number of event triggers before stopping the orbit calculation of
            each particle. Defaults to 0, which makes the event non-terminal.
        """
        logger.info("Calculating particle's orbits...")

        times = []
        reached_end = 0
        terminated = 0

        for p in self.particles:

            logger.remove()

            start = time()
            event = when_psi(p.psi0, terminal)
            p.run(info=False, orbit=orbit, events=[event])
            times.append(time() - start)

            add_logger

            if p.message[0] == "0":
                reached_end += 1
            else:
                terminated += 1

            if "0" == "1":
                pass

        add_logger()

        times = np.array(times)
        time_str = (
            f"Total calculation time: {times.sum():.4g}s.\n"
            + f"Fastest particle: {times.min():.4g}s.\n"
            + f"Slowest particle: {times.max():.4g}s.\n"
            + f"Average time: ({times.mean():.4g} \u00B1 {times.std():.4g})s.\n"
            + f"{reached_end} particles reached the end of the integration interval.\n"
            + f"{terminated} particles were terminated early.\n"
        )
        print(time_str)
        logger.info(time_str)
