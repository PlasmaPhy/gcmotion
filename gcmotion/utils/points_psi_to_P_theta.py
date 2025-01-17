from collections import deque
import numpy as np
from gcmotion.entities.profile import Profile


# Takes in points list/array/deque of the form [theta, psi] returns [theta, P_theat]
def points_psi_to_P_theta(points, profile: Profile):

    new_points = points.copy

    if isinstance(points, np.ndarray):
        new_points = deque(new_points)

    for i in range(len(new_points)):

        theta, psi = new_points[i]

        P_theta = profile.findPtheta(psi=abs(psi))

        new_points[i] = [theta, P_theta]

    return new_points
