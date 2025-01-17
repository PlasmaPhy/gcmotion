import numpy as np
from collections import deque
from gcmotion.entities.profile import Profile
from gcmotion.utils.second_derivative import higher_order_second_derivative
from gcmotion.utils.points_psi_to_P_theta import points_psi_to_P_theta


def XO_points_classification(
    unclassified_fixed_points: np.ndarray,
    profile: Profile,
    delta: float = 1e-5,
    to_P_thetas: bool = True,
):
    # Parameters
    mu = profile.muNU.m
    Pzeta = profile.Pzeta.m

    O_points = deque([])  # Deque for stable O-points
    X_points = deque([])  # Deque for unstable X-points and saddle X-points

    def WNU(theta, psi):
        # Calculate the Hamiltonian at (theta, psi)
        W = profile.findEnergy(psi=abs(psi), theta=theta, units="NUJoule", potential=True)

        return W

    for fixed_point in unclassified_fixed_points:

        theta_fixed, psi_fixed = fixed_point

        # Compute the Hessian matrix elements
        d2W_dtheta2 = higher_order_second_derivative(WNU, theta_fixed, psi_fixed, delta, delta, "x")
        d2W_dPtheta2 = higher_order_second_derivative(
            WNU, theta_fixed, psi_fixed, delta, delta, "y"
        )
        d2W_dtheta_dPtheta = higher_order_second_derivative(
            WNU, theta_fixed, psi_fixed, delta, delta, "mixed"
        )
        # Hessian matrix
        Hessian = np.array([[d2W_dtheta2, d2W_dtheta_dPtheta], [d2W_dtheta_dPtheta, d2W_dPtheta2]])

        # Determinant of the Hessian
        det_Hessian = np.linalg.det(Hessian)

        # Classification based on determinant
        if det_Hessian < 0:
            X_points.append(fixed_point)
        elif det_Hessian > 0:
            O_points.append(fixed_point)

    if to_P_thetas:
        # We have XO points, each, of the form [thetaXO, psiXO] and we will transform them
        # to [thetaXO, P_theta_XO], if asked
        X_points = points_psi_to_P_theta(X_points, profile=profile)
        O_points = points_psi_to_P_theta(O_points, profile=profile)

    return X_points, O_points
