# Import the logger first
from gcmotion.utils import _logger_setup

from gcmotion.classes.particle import Particle
from gcmotion.classes.collection import Collection
from gcmotion.scripts import events

from gcmotion.tokamak import qfactor
from gcmotion.tokamak import bfield
from gcmotion.tokamak import efield


from gcmotion.plotters.time_evolution import time_evolution
from gcmotion.plotters.tokamak_profile import tokamak_profile
from gcmotion.plotters.drift import drift
from gcmotion.plotters.drifts import drifts
from gcmotion.plotters.energy_contour import energy_contour
from gcmotion.plotters.parabolas import parabolas
from gcmotion.plotters.torus2d import torus2d
from gcmotion.plotters.torus3d import torus3d

from gcmotion.plotters.collection_drift import collection_drift
from gcmotion.plotters.collection_energy_contour import collection_energy_contour


import gcmotion.scripts.animation as animation

__all__ = [
    "_logger_setup",
    "Particle",
    "Collection",
    "qfactor",
    "bfield",
    "efield",
    "time_evolution",
    "tokamak_profile",
    "drift",
    "drifts",
    "energy_contour",
    "parabolas",
    "torus2d",
    "torus3d",
    "collection_drift",
    "collection_energy_contour",
    "events",
    "animation",
]
