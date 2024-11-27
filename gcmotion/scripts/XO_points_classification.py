import numpy as np
from collections import namedtuple, deque
from gcmotion.utils.energy_Ptheta import energy_Ptheta
from gcmotion.utils.second_derivative import higher_order_second_derivative
from gcmotion.utils.second_derivative import higher_order_second_derivative


def XO_points_classification(
    unclassified_fixed_points: np.ndarray,
    parameters: namedtuple,
    profile: namedtuple,
    delta: float = 1e-5,
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

    # We have XO points, each, of the form [thetaXO, psiXO] and we will transform them
    # to [thetaXO, P_theta_XO]
    for i in range(len(X_points)):

        X_point = X_points[i]

        theta_X, psi_X = X_point

        _, P_theta_X = energy_Ptheta(
            psi=psi_X, theta=theta_X, mu=mu, Pzeta=Pzeta0, profile=profile, contour_Phi=True
        )

        X_points[i] = [theta_X, P_theta_X]

    for i in range(len(O_points)):

        O_point = O_points[i]

        theta_O, psi_O = O_point

        _, P_theta_O = energy_Ptheta(
            psi=psi_O, theta=theta_O, mu=mu, Pzeta=Pzeta0, profile=profile, contour_Phi=True
        )

        O_points[i] = [theta_O, P_theta_O]

    return X_points, O_points
