from .radial_derivatives_dynamic import (radial_derivatives_liquid_general as radial_derivatives_liquid_dynamic,
                                         radial_derivatives_solid_general as radial_derivatives_solid_dynamic)
from .radial_derivatives_static import (radial_derivatives_liquid_general as radial_derivatives_liquid_static,
                                        radial_derivatives_solid_general as radial_derivatives_solid_static)

from .odes import dynamic_liquid_ode, dynamic_solid_ode, static_liquid_ode, static_solid_ode


# Stored by is_solid, is_static
known_multilayer_odes = {
    (True, True)   : static_solid_ode,
    (True, False)  : dynamic_solid_ode,
    (False, True)  : static_liquid_ode,
    (False, False) : dynamic_liquid_ode
    }