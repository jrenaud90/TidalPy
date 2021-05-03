# Use the Kamata+2015 solution as the default for both liquid and solid dynamic tides. Use Saito for liquid static tides
from .initial_solution_dynamic import solid_guess_kamata as solid_dynamic_guess, \
    liquid_guess_kamata as liquid_dynamic_guess, SolidDynamicGuess, LiquidDynamicGuess
from .initial_solution_static import solid_guess_kamata as solid_static_guess, \
    liquid_guess_saito as liquid_static_guess, SolidStaticGuess, LiquidStaticGuess
from .radial_derivatives_dynamic import radial_derivatives_solid_general as radial_derivatives_solid_dynamic, \
    radial_derivatives_liquid_general as radial_derivatives_liquid_dynamic, RadialFuncSolidDynamicType, \
    RadialFuncLiquidDynamicType
from .radial_derivatives_static import radial_derivatives_solid_general as radial_derivatives_solid_static, \
    radial_derivatives_liquid_general as radial_derivatives_liquid_static, RadialFuncSolidStaticType, \
    RadialFuncLiquidStaticType
