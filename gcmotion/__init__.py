"""
==============================
Guiding Center Motion
==============================


Classes
==============================

.. autosummary::
   :toctree: 

   Particle

Tokamak Configuration
==============================

.. autosummary::
    :toctree: 
    
    qfactor
    bfield
    efield

Plotters
==============================

.. autosummary::
    :toctree: 

    tokamak_profile
    time_evolution
    drift
    drifts
    contour_energy
    parabolas
    torus2d
    torus3d

"""

# Import the logger first
from gcmotion.utils import _logger_setup

from gcmotion.classes.particle import Particle

from gcmotion.tokamak import qfactor
from gcmotion.tokamak import bfield
from gcmotion.tokamak import efield


from gcmotion.plotters.time_evolution import time_evolution
from gcmotion.plotters.tokamak_profile import tokamak_profile
from gcmotion.plotters.drift import drift
from gcmotion.plotters.drifts import drifts
from gcmotion.plotters.contour_energy import contour_energy
from gcmotion.plotters.parabolas import parabolas
from gcmotion.plotters.torus2d import torus2d
from gcmotion.plotters.torus3d import torus3d

__all__ = [
    "_logger_setup",
    "Particle",
    "collection",
    "qfactor",
    "bfield",
    "efield",
    "time_evolution",
    "tokamak_profile",
    "drift",
    "drifts",
    "contour_energy",
    "parabolas",
    "torus2d",
    "torus3d",
]
