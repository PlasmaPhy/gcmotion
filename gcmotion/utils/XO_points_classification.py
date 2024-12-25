import numpy as np
from collections import namedtuple, deque
from gcmotion.utils.energy_Ptheta import energy_Ptheta
from gcmotion.utils.second_derivative import higher_order_second_derivative
from gcmotion.utils.points_psi_to_P_theta import points_psi_to_P_theta


def XO_points_classification(
    unclassified_fixed_points: np.ndarray,
    parameters: namedtuple,
    profile: namedtuple,
    delta: float = 1e-5,
    to_P_thetas: bool = True,
):
    # Parameters
    mu = parameters.mu
    Pzeta0 = parameters.Pzeta0

    O_points = deque([])  # Deque for stable O-points
    X_points = deque([])  # Deque for unstable X-points and saddle X-points

    def WNU(theta, psi):
        # Calculate the Hamiltonian at (theta, psi)
        W, _ = energy_Ptheta(
            psi=psi, theta=theta, Pzeta=Pzeta0, mu=mu, profile=profile, contour_Phi=True
        )

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
        X_points = points_psi_to_P_theta(X_points, Pzeta=Pzeta0, mu=mu, profile=profile)
        O_points = points_psi_to_P_theta(O_points, Pzeta=Pzeta0, mu=mu, profile=profile)

    return X_points, O_points
