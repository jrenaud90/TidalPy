from typing import Union

from TidalPy.utilities.performance import njit

from .functions import takeuchi_phi_psi, takeuchi_phi_psi_general, z_calc

from .initial_solution_dynamic import LiquidDynamicGuess, SolidDynamicGuess
from .initial_solution_static import LiquidStaticGuess, SolidStaticGuess

InitialGuess = Union[LiquidDynamicGuess, LiquidStaticGuess, SolidDynamicGuess, SolidStaticGuess]

from .initial_solution_dynamic import (
    liquid_guess_kamata   as liquid_dynamic_compressible_kmn15,
    liquid_guess_takeuchi as liquid_dynamic_compressible_other,
    solid_guess_kamata    as solid_dynamic_compressible_kmn15,
    solid_guess_takeuchi  as solid_dynamic_compressible_other)

from .initial_solution_dynamic_incomp import (
    liquid_guess_kamata   as liquid_dynamic_incompressible_kmn15,
    liquid_guess_takeuchi as liquid_dynamic_incompressible_other,
    solid_guess_kamata    as solid_dynamic_incompressible_kmn15,
    solid_guess_takeuchi  as solid_dynamic_incompressible_other)

from .initial_solution_static import (
    liquid_guess_saito   as liquid_static_compressible_other,
    solid_guess_kamata   as solid_static_compressible_kmn15,
    solid_guess_takeuchi as solid_static_compressible_other)

from .initial_solution_static_incomp import (
    liquid_guess_saito   as liquid_static_incompressible_other,
    solid_guess_kamata   as solid_static_incompressible_kmn15,
    solid_guess_takeuchi as solid_static_incompressible_other)

# Kamata et al. (2015) did not provide equations for static liquid layers so use the Saito (1974) method instead.
liquid_static_compressible_kmn15   = liquid_static_compressible_other
liquid_static_incompressible_kmn15 = liquid_static_incompressible_other

# Known initial guess functions stored by: is_solid, is_static, is_incompressible, is_kamata
known_initial_guess_funcs = {
    (True, True, True, True)    : solid_static_incompressible_kmn15,
    (True, True, True, False)   : solid_static_incompressible_other,
    (True, True, False, True)   : solid_static_compressible_kmn15,
    (True, True, False, False)  : solid_static_compressible_other,
    (True, False, True, True)   : solid_dynamic_incompressible_kmn15,
    (True, False, True, False)  : solid_dynamic_incompressible_other,
    (True, False, False, True)  : solid_dynamic_compressible_kmn15,
    (True, False, False, False) : solid_dynamic_compressible_other,
    (False, True, True, True)   : liquid_static_incompressible_kmn15,
    (False, True, True, False)  : liquid_static_incompressible_other,
    (False, True, False, True)  : liquid_static_compressible_kmn15,
    (False, True, False, False) : liquid_static_compressible_other,
    (False, False, True, True)  : liquid_dynamic_incompressible_kmn15,
    (False, False, True, False) : liquid_dynamic_incompressible_other,
    (False, False, False, True) : liquid_dynamic_compressible_kmn15,
    (False, False, False, False): liquid_dynamic_compressible_other,
    }


@njit(cacheable=True)
def find_initial_guess(
        is_solid: bool, is_static: bool, is_incompressible: bool, is_kamata: bool,
        radius: float, shear_modulus: Union[float, complex], bulk_modulus: Union[float, complex],
        density: float, frequency: float, order_l: int, G_to_use: float
    ) -> InitialGuess:
    """

    Parameters
    ----------
    is_solid : bool
        Flag for if the layer is solid (True) or liquid (False).
    is_static : bool
        Flag for if the layer is static (True) or dynamic (False).
    is_incompressible : bool
        Flag for if the layer is incompressible (True) or compressible (False).
    is_kamata : bool
        Flag for if initial conditions should be calculated using the equations from Kamata et al. (2015) or other refs.
    radius : float
        Radius where initial guess is calculated [m].
        This radius should not be zero.
    shear_modulus : Union[float, complex]
        Shear modulus (complex or real) at radius [Pa].
    bulk_modulus : Union[float, complex]
        Bulk modulus (complex or real) at radius [Pa].
    density : float
        Density at radius [kg m-3].
    frequency : float
        Tidal forcing frequency [rad s-1].
    order_l : int
        Tidal harmonic degree.
    G_to_use : float
        Gravitational constant. Use non-dimensional version if using calculating non-dimensionalized radial solutions.

    Returns
    -------
    initial_guess : InitialGuess
        Initial (starting) guess at the innermost region of interest (often the center of a planet).
        The initial guess will be a list of solutions (1 to 3 depending on the assumptions).
        Each solution will be a np.ndarray (2 to 6 values depending on the assumptions).
    """

    if is_solid:
        if is_static:
            if is_incompressible:
                if is_kamata:
                    initial_guess = solid_static_incompressible_kmn15(
                            radius, shear_modulus, density,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = solid_static_incompressible_other(
                            radius, shear_modulus, density,
                            order_l=order_l, G_to_use=G_to_use
                            )
            else:
                if is_kamata:
                    initial_guess = solid_static_compressible_kmn15(
                            radius, shear_modulus, bulk_modulus, density,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = solid_static_compressible_other(
                            radius, shear_modulus, bulk_modulus, density,
                            order_l=order_l, G_to_use=G_to_use
                            )
        else:
            if is_incompressible:
                if is_kamata:
                    initial_guess = solid_dynamic_incompressible_kmn15(
                            radius, shear_modulus, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = solid_dynamic_incompressible_other(
                            radius, shear_modulus, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
            else:
                if is_kamata:
                    initial_guess = solid_dynamic_compressible_kmn15(
                            radius, shear_modulus, bulk_modulus, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = solid_dynamic_compressible_other(
                            radius, shear_modulus, bulk_modulus, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
    else:
        if is_static:
            if is_incompressible:
                if is_kamata:
                    initial_guess = liquid_static_incompressible_kmn15(
                            radius,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = liquid_static_incompressible_other(
                            radius,
                            order_l=order_l, G_to_use=G_to_use
                            )
            else:
                if is_kamata:
                    initial_guess = liquid_static_compressible_kmn15(
                            radius,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = liquid_static_compressible_other(
                            radius,
                            order_l=order_l, G_to_use=G_to_use
                            )
        else:
            if is_incompressible:
                if is_kamata:
                    initial_guess = liquid_dynamic_incompressible_kmn15(
                            radius, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = liquid_dynamic_incompressible_other(
                            radius, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
            else:
                if is_kamata:
                    initial_guess = liquid_dynamic_compressible_kmn15(
                            radius, bulk_modulus, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )
                else:
                    initial_guess = liquid_dynamic_compressible_other(
                            radius, bulk_modulus, density, frequency,
                            order_l=order_l, G_to_use=G_to_use
                            )

    return initial_guess
