from collections import namedtuple, deque
import numpy as np
from gcmotion.utils.energy_Ptheta import energy_Ptheta


# Takes in points list/array/deque of the form [theta, psi] returns [theta, P_theat]
def points_to_P_theta(points, Pzeta: float, mu: float, profile: namedtuple):

    if isinstance(points, np.ndarray):
        points = deque(points)

    for i in range(len(points)):

        theta, psi = points[i]

        _, P_theta = energy_Ptheta(
            theta=theta, psi=psi, mu=mu, Pzeta=Pzeta, profile=profile, contour_Phi=True
        )

        points[i] = [theta, P_theta]

    return points
