# TidalPy Major Change Log


## Version 0.2.0 (September 2019)
* Integration module added
    * Time integration studies can now be performed in a Functional-programming scheme (OOP not directly supported).
    * Currently no differential equations are provided, only the tools to solve them. The user must build their own diff-eqs in a separate script and then import the solvers from TidalPy.integration.
    * OOP classes can be used in the solvers but must be wrapped by njit-capable functions.
        * For example: Using a planet core's mass that was calculated via Burnman can be used like:
~~~~
from TidalPy.planets import build_planet
from TidalPy.integration.multi_rheo_integration import multi_rheo_integration
from numba import njit

planet = build_planet('charon')
core_mass = planet.core.mass

@njit
def my_custom_diffeq(time, dependent_variables):
    # Do stuff
    core_density = core_mass / core_volume
    # Do stuff with core_mass (referenced as a global variable)
    return change_in_dependent_variables	
~~~~

## Version 0.1.0 (July 2019)
Main Release - Changes not tracked