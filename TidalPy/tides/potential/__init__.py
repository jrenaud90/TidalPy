from typing import Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray

TidalPotentialOutput = Tuple['FloatArray', 'FloatArray', 'FloatArray', 'FloatArray', 'FloatArray', 'FloatArray']
PotentialTupleModeOutput = Dict[str, Tuple['FloatArray', 'FloatArray', 'FloatArray', 'FloatArray', 'FloatArray', 'FloatArray']]
TidalPotentialModeOutput = Tuple[Dict[str, 'FloatArray'], Dict[str, 'FloatArray'], PotentialTupleModeOutput]

# Minimum difference between spin and orbital frequency before it is treated as zero.
MIN_SPIN_ORBITAL_DIFF = 1.0e-10

from .synchronous_low_e import tidal_potential as tidal_potential_simple
from .nsr_modes_med_eccen_no_obliquity import tidal_potential as tidal_potential_nsr_modes
from .nsr_med_eccen_no_obliquity import tidal_potential as tidal_potential_nsr
from .nsr_med_eccen_med_obliquity import tidal_potential as tidal_potential_obliquity_nsr
from .nsr_med_eccen_gen_obliquity import tidal_potential as tidal_potential_gen_obliquity_nsr
from .nsr_modes_med_eccen_med_obliquity import tidal_potential as tidal_potential_obliquity_nsr_modes
from .nsr_modes_med_eccen_gen_obliquity import tidal_potential as tidal_potential_gen_obliquity_nsr_modes
from .nsr_modes_low_eccen_gen_obliquity import tidal_potential as tidal_potential_gen_obliquity_low_e_nsr_modes