import numpy as np
import pint
from itertools import product

from gcmotion.entities.profile import Profile
from collections import namedtuple
from gcmotion.utils.distinctify import distinctify
from gcmotion.utils.fp_ic_scan import fp_ic_scan as ic_scanner
from gcmotion.utils.single_fixed_point import fixed_point

# Quantity alias for type annotations
type Quantity = pint.UnitRegistry.Quantity


def _set_up_fixed_points(
    profile: Profile,
    method: str,
    theta_lim: list,
    psi_lim: list,
    fp_ic_scan_tol: float,
    ic_theta_grid_density: int,
    ic_psi_grid_density: int,
    random_init_cond: bool = False,
    LAR_thetas: bool = False,
):

    # CAUTION: Here psi lim has already been denormalized from psi_wall and converted to
    # normalized units. This happened in fixed_points()
    psi_min = psi_lim[0]
    psi_max = psi_lim[1]

    theta_min = theta_lim[0]
    theta_max = theta_lim[1]

    known_theta_values = []

    if LAR_thetas and not random_init_cond:

        bounds = [(0.99 * psi_min, 1.01 * psi_max)]

        initial_conditions = np.linspace(psi_min, psi_max, ic_psi_grid_density)
        known_theta_values = np.linspace(-np.pi, np.pi, 3)

        empty_fixed_points = np.empty((len(known_theta_values) * len(initial_conditions), 2))
        empty_fixed_points[:] = np.nan

    elif random_init_cond and not LAR_thetas:

        bounds = [(theta_min, theta_max), (0.99 * psi_min, 1.01 * psi_max)]

        theta_ic = np.linspace(theta_min, theta_max, ic_theta_grid_density)
        psi_ic = np.linspace(psi_min, psi_max, ic_psi_grid_density)

        initial_conditions = [
            [theta_init, psi_init] for theta_init, psi_init in product(theta_ic, psi_ic)
        ]

        empty_fixed_points = np.empty((len(initial_conditions), len(initial_conditions[0])))
        empty_fixed_points[:] = np.nan

    elif not LAR_thetas and not random_init_cond:

        bounds = [(theta_min, theta_max), (0.99 * psi_min, 1.01 * psi_max)]

        initial_conditions = ic_scanner(
            profile=profile,
            method=method,
            theta_grid_density=ic_theta_grid_density,
            psi_grid_density=ic_psi_grid_density,
            psi_lim=psi_lim,
            theta_lim=theta_lim,
            tol=fp_ic_scan_tol,
        )

        # Make sure there are initial conditions. If not, use randoom
        if not initial_conditions:
            print("\n\nDID NOT FIND INITIAL CONDITIONS. USING RANDOM\n\n")

            theta_ic = np.linspace(theta_min, theta_max, 3)
            psi_ic = np.linspace(psi_min, psi_max, 70)

            initial_conditions = [
                [theta_init, psi_init] for theta_init, psi_init in product(theta_ic, psi_ic)
            ]

        empty_fixed_points = np.empty((len(initial_conditions), len(initial_conditions[0])))
        empty_fixed_points[:] = np.nan

    elif LAR_thetas and random_init_cond:
        print("LAR_thetas and random_init_cond can not be true at the same time.")
        return

    return bounds, initial_conditions, empty_fixed_points, known_theta_values


def fixed_points(
    profile: Profile,
    Q: Quantity,
    method: str = "differential evolution",
    theta_lim: list = [-1.01 * np.pi, 1.01 * np.pi],
    psi_lim: list = [0.01, 1.3],
    dist_tol: float = 1e-3,
    fp_ic_scan_tol: float = 5 * 1e-8,
    ic_theta_grid_density: int = 1000,
    ic_psi_grid_density: int = 1000,
    random_init_cond: bool = False,
    info: bool = False,
    ic_info: bool = False,
    LAR_thetas: bool = False,
):

    psi_lim = np.array(psi_lim) * Q("NUpsi_wall")
    psi_lim = psi_lim.to("NUmagnetic_flux").m

    bounds, initial_conditions, fixed_points, known_theta_values = _set_up_fixed_points(
        profile=profile,
        method=method,
        theta_lim=theta_lim,
        psi_lim=psi_lim,
        fp_ic_scan_tol=fp_ic_scan_tol,
        ic_theta_grid_density=ic_theta_grid_density,
        ic_psi_grid_density=ic_psi_grid_density,
        random_init_cond=random_init_cond,
        LAR_thetas=LAR_thetas,
    )

    idx = 0

    # Run fixed_point() for multiple initial conditions in order to locate
    # multiple fixed points. If thetas=LAR thetas run for the multiple LAR thetas.
    if len(known_theta_values) == 0:
        for initial_condition in initial_conditions:

            theta_fix, psi_fix = fixed_point(
                method=method,
                initial_condition=initial_condition,
                bounds=bounds,
                profile=profile,
                known_thetas=LAR_thetas,
            )
            fixed_points[idx] = [float(theta_fix), float(psi_fix)]
            idx += 1

    else:
        for known_theta_value in known_theta_values:
            for initial_condition in initial_conditions:
                theta_fix, psi_fix = fixed_point(
                    method=method,
                    initial_condition=initial_condition,
                    bounds=bounds,
                    profile=profile,
                    known_thetas=LAR_thetas,
                    known_theta_value=known_theta_value,
                )

                fixed_points[idx] = [float(theta_fix), float(psi_fix)]
                idx += 1

    # A lot of the fixed points that were found have identical values-->
    # find out how many distinct fixed points were located
    distinct_fixed_points = distinctify(fixed_points, tol=dist_tol)
    num_of_dfp = distinct_fixed_points.shape[0]

    if ic_info and not LAR_thetas:
        initial_conditions_print = [[float(x), float(y)] for x, y in initial_conditions]
        print(f"\nInitial Conditions:{initial_conditions_print}\n")
        print(f"\nNumber of Initial Conditions: {len(initial_conditions_print)}\n")

    if info and fixed_points.shape[0] <= 30:
        print(f"\nFixed Points: {fixed_points}\n")
        print(f"Number of Fixed Points: {fixed_points.shape[0]}")

    if info and distinct_fixed_points.shape[0] <= 12:
        print(f"\nDistinct Fixed Points: {distinct_fixed_points}\n")
        print(f"Number of Distinct Fixed Points: {distinct_fixed_points.shape[0]}\n")

    return num_of_dfp, distinct_fixed_points, initial_conditions
