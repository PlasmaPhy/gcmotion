r"""
Simple function wrapper around SciPy's solve_ivp() solver. It solves
the dynamical system:

.. todo::
    link to system formula

It uses the RK4(5) solving method. The tolerance can be tweaked 
in gcmotion/configuration/solver_configuration.py.

The solver only calculates the time evolution of
:math:`\theta, \psi, \zeta, \rho_{||}`. The rest of the dynamical variables
are calculated afterwards.

.. tip::
    When setting up the analytical function dSdt to pass to the solver,
    it is better to use Python's built-in mathematical functions by importing
    the ``math`` module, instead of their ``numpy``
    counterparts. Python's functions are much faster than numpy's when calculating
    single values (often by a factor of 10 or more).

Example
-------

This is how ``orbit`` is called inside the class :py:class:`Particle`:

.. code-block:: python
    
    t = self.t_eval

    init_cond = {
        "theta0": self.theta0,
        "psi0": self.psi0,
        "zeta0": self.zeta0,
        "rho0": self.rho0,
    }

    constants = {
        "mu": self.mu,
    }

    profile = {
        "qfactor": self.qfactor,
        "Bfield": self.Bfield,
        "Efield": self.Efield,
        "Volts_to_NU": self.Volts_to_NU,
    }

    solution orbit(t, init_cond, constants, profile, events)

And this is how the solution is unpacked:

.. code-block:: python

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


.. rubric:: Function:
    :heading-level: 4
"""

import numpy as np
from scipy.integrate import solve_ivp
from math import sqrt

from gcmotion.utils._logger_setup import logger

from gcmotion.configuration.solver_configuration import (
    solver_configuration as config,
)


def orbit(
    t: np.ndarray,
    init_cond: list,
    constants: dict,
    profile: dict,
    events: list,
):
    r"""Wrapper function around SciPy's solve_ivp().

    Parameters
    ----------

    t : np.ndarray
        The evaluation times array.
    init_cond : list
        List containing the initial conditions
        :math:`\theta_0, \psi_0, \zeta_0, \rho_{||,0}`.
    constants : dict
        Dict containing the constants of motion. Currently only :math:`\mu`
        is used.
    profile : dict
        Dict containing the tokamak configuration objects.
    events : list
        List containing the independed events to track.

    Returns
    -------

    solution : dict
        A dict containing the following:

            #. The calculated solutions of
                :math:`\theta, \psi, \zeta, \rho, \psi_p, P_\theta, P_\zeta`
                as np.ndarrays.

            #. A np.ndarray containing the evaluation times. This is useful in the case
                a terminating event is used.

            #. 2 lists containing the time positions and values of the found events. If
                more than 1 event is used, the lists are nested.

            #. The solver status message.

    """
    logger.info(f"Calculating orbit with events {events}")

    # Tokamak profile
    qfactor = profile["qfactor"]
    Bfield = profile["Bfield"]
    Efield = profile["Efield"]
    Volts_to_NU = float(profile["Volts_to_NU"])

    # Constants of motion
    # E = float(constants["E"]) # Not used
    mu = float(constants["mu"])
    mi = float(constants["mi"])
    qi = int(constants["qi"])
    # Pzeta0 = float(constants["Pzeta0"]) # Not used

    # Initial Conditions
    S0 = [
        float(init_cond["theta0"]),
        float(init_cond["psi0"]),
        float(init_cond["zeta0"]),
        float(init_cond["rho0"]),
    ]

    def dSdt(t, S):
        """Sets the diff equations system to pass to scipy.

        All values are in normalized units (NU).
        """

        theta, psi, z, rho = S

        # Intermediate values
        phi_der_psip, phi_der_theta = Efield.Phi_der(psi)
        phi_der_psip *= Volts_to_NU
        phi_der_theta *= Volts_to_NU
        B_der_psi, B_der_theta = Bfield.B_der(psi, theta)
        q = qfactor.q_of_psi(psi)
        r = sqrt(2 * psi)
        B = Bfield.B(r, theta)
        par = mu + (qi**2 / mi) * rho**2 * B
        bracket1 = -par * q * B_der_psi + qi * phi_der_psip
        bracket2 = par * B_der_theta + qi * phi_der_theta
        D = Bfield.g * q + Bfield.I

        # Canonical Equations
        theta_dot = (qi / mi) / D * rho * B**2 + Bfield.g / (D * qi) * bracket1
        psi_dot = -Bfield.g * q / (D * qi) * bracket2
        rho_dot = psi_dot / (Bfield.g * q)
        z_dot = q * (
            (qi / mi) * rho * B**2 / D - Bfield.I / (D * qi) * bracket1
        )

        return [theta_dot, psi_dot, z_dot, rho_dot]

    t_span = (t[0], t[-1])
    sol = solve_ivp(
        dSdt,
        t_span=t_span,
        y0=S0,
        t_eval=t,
        atol=config["atol"],
        rtol=config["rtol"],
        events=events,
        dense_output=True,
    )

    theta = sol.y[0]
    psi = sol.y[1]
    zeta = sol.y[2]
    rho = sol.y[3]
    t_eval = sol.t
    t_events = sol.t_events
    y_events = sol.y_events
    message = f"{sol.status}: {sol.message}"

    logger.debug(f"Solver status: {message}")

    # Calculate psip and Canonical Momenta
    psip = qfactor.psip_of_psi(psi)
    Ptheta = psi + rho * Bfield.I
    Pzeta = rho * Bfield.g - psip

    solution = {
        "theta": theta,
        "psi": psi,
        "zeta": zeta,
        "rho": rho,
        "psip": psip,
        "Ptheta": Ptheta,
        "Pzeta": Pzeta,
        "t_eval": t_eval,
        "t_events": t_events,
        "y_events": y_events,
        "message": message,
    }

    return solution
