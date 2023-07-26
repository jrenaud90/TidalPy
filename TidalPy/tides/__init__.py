from .love1d import (complex_love as calc_complex_love, complex_love_general as calc_complex_love_general,
                     effective_rigidity as calc_effective_rigidity,
                     effective_rigidity_general as calc_effective_rigidity_general,
                     static_love as calc_static_love, static_love_general as calc_static_love_general)
from .methods import GlobalApproxTides, LayeredTides, TidesBase
from .dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced
from .heating import calculate_volumetric_heating