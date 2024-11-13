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
            psi_lim=psi_lim,
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
    psi_lim : list, optional
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
from collections import deque, namedtuple

from gcmotion.scripts.fixed_points import fixed_points
from gcmotion.classes.collection import Collection


def bifurcation(
    collection: Collection,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: list = [0.01, 1.3],
    iterations: int = 250,
    info: bool = False,
):

    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

    num_of_fp = deque([])
    fp = deque([])

    thetas_fixed = deque([])
    P_thetas_fixed = deque([])

    p1 = collection[0]
    p_last = collection[-1]
    Q = p1.Q

    # Check if the partcles have different Pzeta0's
    if p1.Pzeta0 == p_last.Pzeta0:
        print(r"Each particle in the collection must have different $P_{\zeta0}$")
        return

    # Get Tokamak profile
    qfactor = p1.qfactor
    bfield = p1.bfield
    efield = p1.efield

    Profile = namedtuple("Tokamak_Profile", ["qfactor", "bfield", "efield"])
    profile = Profile(
        qfactor=qfactor,
        bfield=bfield,
        efield=efield,
    )

    Parameters = namedtuple("Orbit_Parameters", ["Pzeta0", "mu"])

    for idx, p in enumerate(collection):

        current_P_zeta = p.Pzeta0NU.magnitude
        # Get Particle Parameters
        parameters = Parameters(
            Pzeta0=current_P_zeta,
            mu=p1.muNU.magnitude,
        )

        current_num_of_fp, current_fp = fixed_points(
            parameters=parameters,
            profile=profile,
            Q=Q,
            theta_lim=[theta_min, theta_max],
            psi_lim=psi_lim,
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
