# TidalPy Major Change Log


## Version 0.2.0 (Fall 2019)
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

* Modified OOP Backend
    * Models are now free to access the state properties of layers and planets. This eliminates the need to pass arguments into model classes for calculations. All the user has to do is change the layer and/or planet's state properties and then those changes will automatically propagate.

* Redesigned Dynamics Module
    * All tidal mode information is stored here.
    * The big change in this version is that higher order tidal order `l`, as well as orbital, truncations have been implemented.

* New Tides Module
    * This module contains functionality to calculate tidal heating and torque based on a complex compliance function and the various thermal and orbital parameters.
    * The functions to calculate the Love number are also stored here.

* New Rheology Scheme
    * To better follow real physics, all strength-based models have been moved under the redesigned `Rheology` class.
    * This includes: liquid and solid viscosity calculations, partial melting, and complex compliance
    * Love numbers are also calculated by the `Rheology` class but their functionality is now stored in the `tides` module.

* QOL Improvements
    * Many new docstrings, type hints, and overall clean up of functions and classes. All docstrings now follow the numpy format.
    * Many new tests for both the functional and OOP versions of TidalPy.
    * More comments and spelling/typo fixes everywhere.
        
* Other
    * Various functions and classes have been reworked or removed to simplify the code base. 

## Version 0.1.0 (July 2019)
Main Release - Changes not tracked