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
            dist_tol=dist_tol,
            ic_theta_grid_dencity = 800,
            ic_psi_grid_density = 1200,
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
        to the :math:`\psi` variable. It will be passed into the "bounds" argument of
        :py:func:`fixed_points`. 
    dist_tol : float, optional
        Tolerance that determines distinct fixed points. If both :math:`\theta` and
        :math:`\psi` elements of a fixed point are less than :py:data:`dist_tol` apart
        the two fixed points are not considered distinct.
    ic_theta_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\theta` variable 
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution` 
        will be conducted. Will be passed to :py:func:`fixed_points`.
    ic_psi_grid_density : int, optional
        Integer dictating the theta density with regard to the :math:`\psi` variable 
        of the grid upon which the search for initial conditions for the :py:func:`differential_evolution` 
        will be conducted.  Will be passed to :py:func:`fixed_points`.
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

from gcmotion.utils.XO_points_classification import XO_points_classification as xoc
from gcmotion.scripts.fixed_points import fixed_points
from gcmotion.classes.collection import Collection


def bifurcation(
    collection: Collection,
    theta_lim: list = [-np.pi, np.pi],
    psi_lim: list = [0.01, 1.3],
    dist_tol: float = 1e-3,
    ic_theta_grid_density: int = 1000,
    ic_psi_grid_density: int = 1000,
    random_fp_init_cond: bool = False,
    info: bool = False,
    fp_ic_info: bool = False,
):

    num_of_XP = deque([])
    num_of_OP = deque([])

    X_points = deque([])
    O_points = deque([])

    X_thetas = deque([])
    X_P_thetas = deque([])

    O_thetas = deque([])
    O_P_thetas = deque([])

    p1 = collection[0]
    p_last = collection[-1]

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
            Q=p.Q,
            theta_lim=theta_lim,
            psi_lim=psi_lim,
            dist_tol=dist_tol,
            ic_theta_grid_density=ic_theta_grid_density,
            ic_psi_grid_density=ic_psi_grid_density,
            random_init_cond=random_fp_init_cond,
            info=info,
            ic_info=fp_ic_info,
        )

        # CAUTION: The xoc function takes in psis_fixed but returns P_thetas_fixed
        current_X_points, current_O_points = xoc(
            unclassified_fixed_points=current_fp,
            parameters=parameters,
            profile=profile,
        )

        # Convert deque to numpy arrays for easy manipulation
        current_X_thetas, current_X_P_thetas = (
            zip(*current_X_points) if current_X_points else ([], [])
        )
        current_O_thetas, current_O_P_thetas = (
            zip(*current_O_points) if current_O_points else ([], [])
        )

        if info:
            print(
                f"\nCurrent Step: {idx+1} at P_z = {current_P_zeta} with {current_num_of_fp} fixed points"
            )
            print(f"Current Fixed Points: {current_fp}\n")
            print(
                f"Current X Points: {[[float(thetaX),float(P_thetaX)] for thetaX,P_thetaX in current_X_points]}\n"
            )
            print(
                f"Current O Points: {[[float(thetaO),float(P_thetaO)] for thetaO,P_thetaO in current_O_points]}\n"
            )

        num_of_XP.append(len(current_X_points))
        num_of_OP.append(len(current_O_points))

        X_points.append(current_X_points)
        O_points.append(current_O_points)

        X_thetas.append(current_X_thetas)
        X_P_thetas.append(current_X_P_thetas)

        O_thetas.append(current_O_thetas)
        O_P_thetas.append(current_O_P_thetas)

    return (
        X_thetas,
        X_P_thetas,
        O_thetas,
        O_P_thetas,
        num_of_XP,
        num_of_OP,
    )
