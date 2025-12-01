from .love1d import complex_love
from .love1d import complex_love_general
from .love1d import effective_rigidity
from .love1d import effective_rigidity_general
from .love1d import static_love
from .love1d import static_love_general

from .methods import TidesBase as TidesBase
from .methods import GlobalApproxTides as GlobalApproxTides
from .methods import LayeredTides as LayeredTides

from .dissipation import calc_tidal_susceptibility as calc_tidal_susceptibility
from .dissipation import calc_tidal_susceptibility_reduced as calc_tidal_susceptibility_reduced

from .heating import calculate_volumetric_heating as calculate_volumetric_heating
# Alias functions

calc_complex_love         = complex_love
calc_complex_love_general = complex_love_general
calc_effective_rigidity   = effective_rigidity
calc_static_love          = static_love
calc_static_love_general  = static_love_general
calc_effective_rigidity_general = effective_rigidity_general
