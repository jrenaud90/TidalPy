from typing import Dict, Tuple

from ...utilities.types import FloatArray

TidalPotentialOutput = Tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]

PotentialTupleModeOutput = Dict[str, Tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]]
TidalPotentialModeOutput = Tuple[Dict[str, FloatArray], Dict[str, FloatArray], PotentialTupleModeOutput]
MIN_SPIN_ORBITAL_DIFF = 1.0e-10

from .synchronous_low_e import tidal_potential as tidal_potential_simple
from .nsr_modes_med_eccen_no_obliquity import tidal_potential as tidal_potential_nsr_modes
from .nsr_med_eccen_no_obliquity import tidal_potential as tidal_potential_nsr
from .nsr_med_eccen_med_obliquity import tidal_potential as tidal_potential_obliquity_nsr
from .nsr_modes_med_eccen_med_obliquity import tidal_potential as tidal_potential_obliquity_nsr_modes
