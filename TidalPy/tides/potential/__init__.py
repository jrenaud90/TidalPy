from typing import Tuple

from ...utilities.types import FloatArray

TidalPotentialOutput = Tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]
MIN_SPIN_ORBITAL_DIFF = 1.0e-10

from .synchronous_low_e import tidal_potential as tidal_potential_simple
from .nsr_obliquity_low_e import tidal_potential as tidal_potential_obliquity_nsr
from .nsr_med_e import tidal_potential as tidal_potential_nsr
from .nsr_med_e_modes import tidal_potential as tidal_potential_nsr_modes
