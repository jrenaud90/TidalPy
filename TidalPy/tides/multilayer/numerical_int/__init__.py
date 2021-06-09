# Use the Kamata+2015 solution as the default for both liquid and solid dynamic tides. Use Saito for liquid static tides
from .initial_solution_dynamic import solid_guess_kamata as solid_dynamic_guess, \
    liquid_guess_kamata as liquid_dynamic_guess,\
    solid_guess_takeuchi as solid_dynamic_guess_ts72,\
    liquid_guess_takeuchi as liquid_dynamic_guess_ts72,\
    SolidDynamicGuess, LiquidDynamicGuess
from .initial_solution_static import solid_guess_kamata as solid_static_guess, \
    liquid_guess_saito as liquid_static_guess,\
    solid_guess_takeuchi as solid_static_guess_ts72,\
    SolidStaticGuess, LiquidStaticGuess
from .radial_derivatives_dynamic import radial_derivatives_solid_general as radial_derivatives_solid_dynamic, \
    radial_derivatives_liquid_general as radial_derivatives_liquid_dynamic
from .radial_derivatives_static import radial_derivatives_solid_general as radial_derivatives_solid_static, \
    radial_derivatives_liquid_general as radial_derivatives_liquid_static

# Known initial guess functions stored by: is_kamata, is_solid, is_dynamic
known_initial_guess_funcs = {
    (True, True, True): solid_dynamic_guess,
    (False, True, True): solid_dynamic_guess_ts72,
    (True, False, True): liquid_dynamic_guess,
    (False, False, True): liquid_dynamic_guess_ts72,
    (True, True, False): solid_static_guess,
    (False, True, False): solid_static_guess_ts72,
    # The next two are identical because they both come from Saito74 not TS72 or Kamata
    (True, False, False): liquid_static_guess,
    (False, False, False): liquid_static_guess
}

def find_initial_guess(is_kamata: bool, is_solid: bool, is_dynamic: bool):
    # Solve the function too? No, for now...
    return known_initial_guess_funcs[(is_kamata, is_solid, is_dynamic)]