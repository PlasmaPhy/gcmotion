from collections import deque
import numpy as np
from gcmotion.entities.profile import Profile


# Takes in points list/array/deque of the form [theta, psi] returns [theta, P_theat]
def points_psi_to_P_theta(points, profile: Profile):

    # Define a Quantity object. Will be needed later
    Q = profile.Q

    new_points = points.copy()

    if isinstance(points, np.ndarray):
        new_points = deque(new_points)

    for i in range(len(new_points)):

        theta, psi = new_points[i]

        P_thetaNU = profile.findPtheta(psi=Q(abs(psi), "NUMagnetic_flux"))

        new_points[i] = [theta, P_thetaNU.m]

    return new_points
