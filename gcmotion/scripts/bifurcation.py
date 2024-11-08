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
