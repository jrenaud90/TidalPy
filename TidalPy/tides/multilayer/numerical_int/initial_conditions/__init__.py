from .functions import takeuchi_phi_psi, takeuchi_phi_psi_general, z_calc

from .initial_solution_dynamic import (liquid_guess_kamata as liquid_dynamic_guess_kmn15,
                                       liquid_guess_takeuchi as liquid_dynamic_guess_ts72,
                                       solid_guess_kamata as solid_dynamic_guess_kmn15,
                                       solid_guess_takeuchi as solid_dynamic_guess_ts72)
from .initial_solution_static import (liquid_guess_saito as liquid_static_guess_s74,
                                      solid_guess_kamata as solid_static_guess_kmn15,
                                      solid_guess_takeuchi as solid_static_guess_ts72)
from .....utilities.performance import njit

# Known initial guess functions stored by: is_kamata, is_solid, is_dynamic
known_initial_guess_funcs = {
    (True, True, True)   : solid_dynamic_guess_kmn15,
    (False, True, True)  : solid_dynamic_guess_ts72,
    (True, False, True)  : liquid_dynamic_guess_kmn15,
    (False, False, True) : liquid_dynamic_guess_ts72,
    (True, True, False)  : solid_static_guess_kmn15,
    (False, True, False) : solid_static_guess_ts72,
    # The next two are identical because they both come from Saito74 not TS72 or Kamata
    (True, False, False) : liquid_static_guess_s74,
    (False, False, False): liquid_static_guess_s74
    }

@njit(cacheable=True)
def find_initial_guess(is_kamata: bool, is_solid: bool, is_dynamic: bool,
                       radius, shear_modulus, bulk_modulus, density, frequency, order_l, G_to_use):

    if is_solid:
        if is_dynamic:
            if is_kamata:
                result = solid_dynamic_guess_kmn15(radius, shear_modulus, bulk_modulus, density, frequency,
                                                   order_l=order_l, G_to_use=G_to_use)
            else:
                result = solid_dynamic_guess_ts72(radius, shear_modulus, bulk_modulus, density, frequency,
                                                  order_l=order_l, G_to_use=G_to_use)
        else:
            if is_kamata:
                result = solid_static_guess_kmn15(radius, shear_modulus, bulk_modulus, density,
                                                  order_l=order_l, G_to_use=G_to_use)
            else:
                result = solid_static_guess_ts72(radius, shear_modulus, bulk_modulus, density,
                                                  order_l=order_l, G_to_use=G_to_use)
    else:
        if is_dynamic:
            if is_kamata:
                result = liquid_dynamic_guess_kmn15(radius, bulk_modulus, density, frequency,
                                                    order_l=order_l, G_to_use=G_to_use)
            else:
                result = liquid_dynamic_guess_ts72(radius, bulk_modulus, density, frequency,
                                                   order_l=order_l, G_to_use=G_to_use)
        else:
            # The next two are identical because they both come from Saito74 not TS72 or Kamata
            if is_kamata:
                result = liquid_static_guess_s74(radius, order_l=order_l, G_to_use=G_to_use)
            else:
                result = liquid_static_guess_s74(radius, order_l=order_l, G_to_use=G_to_use)

    return result