import numpy as np

from TidalPy.constants import G

from .radial_derivatives_dynamic import (radial_derivatives_liquid_general as radial_derivatives_liquid_dynamic,
                                         radial_derivatives_solid_general as radial_derivatives_solid_dynamic)
from .radial_derivatives_static import (radial_derivatives_liquid_general as radial_derivatives_liquid_static,
                                        radial_derivatives_solid_general as radial_derivatives_solid_static)
from .radial_derivatives_dynamic_incomp import (
    radial_derivatives_liquid_general as radial_derivatives_liquid_dynamic_incomp,
    radial_derivatives_solid_general as radial_derivatives_solid_dynamic_incomp)
from .radial_derivatives_static_incomp import (
    radial_derivatives_liquid_general as radial_derivatives_liquid_static_incomp,
    radial_derivatives_solid_general as radial_derivatives_solid_static_incomp)

from .odes import dynamic_liquid_ode, dynamic_solid_ode, static_liquid_ode, static_solid_ode

# Known initial guess functions stored by: is_solid, is_static
# Incompressible flag is built into the outer odes.
known_multilayer_odes = {
    (True, True)  : static_solid_ode,
    (True, False) : dynamic_solid_ode,
    (False, True) : static_liquid_ode,
    (False, False): dynamic_liquid_ode
    }


def find_ode(is_solid: bool, is_static: bool, is_incompressible: bool,
             radius_array: np.ndarray, shear_modulus_array, bulk_modulus_array,
             density_array: np.ndarray, gravity_array: np.ndarray, frequency: float,
             order_l: int = 2, G_to_use: float = G) -> tuple:
    """ Find additional arguments that are required for the ODE functions.

    Parameters
    ----------
    is_solid : bool
        Flag for if the layer is solid (true) or liquid (false).
    is_static : bool
        Flag for if the layer is static (true) or dynamic (false).
    is_incompressible : bool
        If `True`, the incompressible assumption will be used.
    radius_array
        Array of radii for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radii of a planet should not be equal to zero [m].
    shear_modulus_array
        Shear modulus at each `radii` [Pa] (can be complex for shear dissipation).
    bulk_modulus_array
        Bulk modulus at each `radii` [Pa] (can be complex for bulk dissipation).
    density_array
        Density at each `radii` [kg m-3].
    gravity_array
        Acceleration due to gravity calculated at each `radii` [m s-2].
    frequency : float
        Forcing frequency [rad s-1].
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    ode : callable
        Differential equation function based on the provided assumptions.
    additional_arg_tuple : tuple
        Additional inputs that the selected ode requires.

    """
    
    ode = known_multilayer_odes[is_solid, is_static]
    if is_solid:
        if is_static:
            additional_arg_tuple = (radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array,
                                    order_l, G_to_use, is_incompressible)
        else:
            additional_arg_tuple = (radius_array, shear_modulus_array, bulk_modulus_array, density_array, gravity_array,
                                    frequency,
                                    order_l, G_to_use, is_incompressible)
    else:
        if is_static:
            additional_arg_tuple = (radius_array, density_array, gravity_array,
                                    order_l, G_to_use, is_incompressible)
        else:
            additional_arg_tuple = (radius_array, bulk_modulus_array, density_array, gravity_array, frequency,
                                    order_l, G_to_use, is_incompressible)
    
    return ode, additional_arg_tuple
