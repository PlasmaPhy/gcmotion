import numpy as np
import matplotlib.pyplot as plt
from math import isclose

from matplotlib.path import Path

from gcmotion.utils.logger_setup import logger
from gcmotion.scripts.utils.shoelace import shoelace
from gcmotion.configuration.scripts_configuration import ContourFreqConfig

config = ContourFreqConfig()
tau = 2 * np.pi


class ContourSegment:

    def __init__(self, segment: Path, profile: float, psilim: tuple, E):
        self.vertices = segment.vertices
        self.psilim = psilim
        self.E = E

        # Bounding box
        [self.theta_min, self.psi_min], [self.theta_max, self.psi_max] = (
            segment.get_extents().get_points()
        )

        # Run in that order
        self.validate()
        if self.valid:
            self.classify()

        # if self.valid:
        #     self.calculate_area()

    def validate(self) -> None:
        r"""Checks if the bbox of the contour line touches the upper or lower
        walls, which means the path gets cut off and must be discarded."""
        inbounds = _is_inbounds(self)
        is_cutoff_trapped = _is_cutoff_trapped(self)
        self.valid = inbounds and not is_cutoff_trapped

    def classify(self):
        self.passing, self.trapped = _tp_classify(self)
        self.left_to_right = _is_left_to_right(self)

    def calculate_area(self):
        extra_left = [[2 * np.pi, 0], [-2 * np.pi, 0]]
        extra_right = extra_left[::-1]

        if self.passing and self.left_to_right:
            first = self.vertices[0]
            extra = extra_left + [first]
            self.vertices = np.append(self.vertices, extra, axis=0)
        elif self.passing and not self.left_to_right:
            last = self.vertices[0]  # right-to-left!
            extra = extra_right + [last]
            self.vertices = np.append(self.vertices, extra, axis=0)
        elif self.trapped:
            pass
        else:
            return None

        area = shoelace(self.vertices)
        if self.passing:
            area /= 2  # because theta span = 4π
        self.J = area / (2 * np.pi)

        return self.J

    def bbox_distance(self, xy: tuple[float, float]):
        r"""Returns the distance (squared) of the origin point from self's
        origin point."""
        return (self.theta_min - xy[0]) ** 2 + (self.psi_min - xy[1]) ** 2

    def simple_plot(self, ax=None, show=True):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot()

        # color = "b" if self.valid else "r"
        color = "b" if self.passing else "r"

        x, y = self.vertices.T[:]
        ax.plot(x, y, color=color, linewidth=1)
        if show:
            plt.show()


# ================================ Validation ================================


def _is_inbounds(path: ContourSegment) -> bool:
    r"""Checks if the path's bounding box is in bounds of the whole contour,
    e.g the contour line doesn't get cutoff.
    """
    atol = 1e-7  # Must not be 0(default) when comparing with 0
    return not (
        isclose(path.psi_min, path.psilim[0], abs_tol=atol)  # Touches floor
        or isclose(path.psi_max, path.psilim[1], abs_tol=atol)  # Touches ceil
    )


def _is_cutoff_trapped(path: ContourSegment) -> bool:
    r"""Checks that the bounding box doesn't touch the left *and* the right
    wall, but touches at least one of them.
    """
    return (
        isclose(path.theta_min, -tau) and not isclose(path.theta_max, tau)
    ) or (isclose(path.theta_max, tau) and not isclose(path.theta_min, -tau))


# ===================== Trapped - Passing Classification =====================


def _tp_classify(path: ContourSegment) -> list[bool, bool]:
    r"""Checks the left and right bbox edges to check if passing or trapped."""

    if isclose(path.theta_min, -tau) and isclose(path.theta_max, tau):
        return True, False
    else:
        return False, True


# ======================= Left-to-Right Classification =======================


def _is_left_to_right(path: ContourSegment) -> bool:
    r"""In the case that the path is passing, it returns True if the *first*
    point touches the left wall. It also checks that the *last* point touches
    the third wall, which might be reduntant."""
    if path.trapped:
        return None

    return (
        True
        if (
            isclose(path.vertices[0][0], -tau)
            and isclose(path.vertices[-1][0], tau)
        )
        else False
    )
