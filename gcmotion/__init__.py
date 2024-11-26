# Import the logger first
from gcmotion.utils import logger_setup

from gcmotion.utils.setup_pint import setup_pint
from gcmotion.utils.get_size import get_size

# main namespace
from gcmotion.tokamak import qfactor
from gcmotion.tokamak import bfield
from gcmotion.tokamak import efield

from gcmotion.entities.tokamak import Tokamak
from gcmotion.entities.physical_parameters import PhysicalParameters
from gcmotion.entities.initial_conditions import InitialConditions

from gcmotion.entities.profile import Profile
from gcmotion.entities.particle import Particle


# main namespace
__all__ = [
    "setup_pint",
    "get_size",
    # Tokamak objects
    "qfactor",
    "bfield",
    "efield",
    # Base entities
    "Tokamak",
    "PhysicalParameters",
    "InitialConditions",
    # Entities
    "Profile",
    "Particle",
]
