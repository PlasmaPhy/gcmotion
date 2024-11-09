r""" Function that calculates the fixed points and the number of fixed points for multiple
particles in a Collection form :py:class:`~gcmotion.classes.collection.Collection`, where 
each particle has a different :math:`P_{\zeta0}`.

Example
-------

This is how :py:func:`bifurcation` can be called inside the function :py:func:`bifurcation_plot`:

.. code-block:: python

    from gcmotion.scripts.bifurcation import bifurcation

        thetas_fixed, P_thetas_fixed, num_of_fp = bifurcation(
            collection=collection,
            theta_lim=theta_lim,
            P_theta_lim=P_theta_lim,
            iterations=iterations,
            info=info,
        )

    Parameters
    ----------
    
    collection : :py:class:`~gcmotion.classes.collection.Collection`
        The collection of particles
    theta_lim : list, optional
        Provides the limits for the solution search area for fixed points
        with regards to the :math:`\theta` variable. It will be passed into 
        :py:func:`fixed_points` 
    P_theta_lim : list, optional
        Provides the limits (divided by psi_wall) for the solution search area with regards
        to the :math:`P_{\theta}` variable. It will be passed into the "bounds" argument of
        :py:func:`fixed_points`. 
    iterations : int, optional
        Integer that essentially dictates the number of initial conditions that will be 
        passed into :py:func:`fixed_points`.
    info : bool, optional
        Boolean that dictates weather the :math:`P_{\zeta0}` of the particle whose
        fixed points have just been calculated, will be printed alongside the fixed points
        found.

    Returns
    -------

    thetas_fixed, P_thetas_fixed, num_of_fp : tuple
        Tuple where each element is a list containing the lists of all the :math:`theta`'s  
        fixed, all the :math:`P_{theta}`'s fixed and the number of fixed points found for 
        each :math:`P_{\zeta0}`.
"""

import numpy as np
from collections import deque

from gcmotion.scripts.fixed_points import fixed_points
from gcmotion.classes.collection import Collection


def bifurcation(
    collection: Collection,
    theta_lim: list = [-np.pi, np.pi],
    P_theta_lim: list = [0.01, 1.3],
    iterations: int = 250,
    info: bool = False,
):

    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

    P_theta_min = P_theta_lim[0]
    P_theta_max = P_theta_lim[1]

    num_of_fp = deque([])
    fp = deque([])

    thetas_fixed = deque([])
    P_thetas_fixed = deque([])

    particles = collection.particles
    p1 = particles[0]
    p_last = particles[-1]

    if p1.Pzeta0 == p_last.Pzeta0:
        print(r"Each particle in the collection must have different $P_{\zeta0}$")
        return

    qfactor = p1.qfactor
    Bfield = p1.Bfield
    Efield = p1.Efield
    Volts_to_NU = float(p1.Volts_to_NU)

    profile = {"qfactor": qfactor, "Bfield": Bfield, "Efield": Efield, "Volts_to_NU": Volts_to_NU}

    mu = float(p1.mu)
    mi = float(p1.mi)
    qi = float(p1.qi)
    Pzeta0 = float(p1.Pzeta0)

    constants = {"mu": mu, "mass": mi, "qi": qi, "Pzeta0": Pzeta0}

    for idx, p in enumerate(particles):

        current_P_zeta = p.Pzeta0
        constants["Pzeta0"] = current_P_zeta
        psi_wall = p.psi_wall

        current_num_of_fp, current_fp = fixed_points(
            profile=profile,
            constants=constants,
            theta_lim=[theta_min, theta_max],
            P_theta_lim=[P_theta_min * psi_wall, P_theta_max * psi_wall],
            iterations=iterations,
            info=False,
        )

        current_thetas_fixed = current_fp[:, 0]
        current_P_thetas_fixed = current_fp[:, 1]

        if info:
            print(
                f"\nCurrent Step: {idx+1} at P_z = {current_P_zeta} with {current_num_of_fp} fixed points"
            )
            print(f"Current Fixed Points: {current_fp}\n")

        fp.append(current_fp)
        num_of_fp.append(current_num_of_fp)

        thetas_fixed.append(current_thetas_fixed)
        P_thetas_fixed.append(current_P_thetas_fixed)

    return thetas_fixed, P_thetas_fixed, num_of_fp
