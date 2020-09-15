# TidalPy Major Change Log

## Version 0.2.1 Alpha (Fall 2019)
*Will break studies based on previous versions*

Note: TidalPy version of "0.2.0" was never made publicly available. This is the next version after 0.1.0.

* Integration module added
    * Time integration studies can now be performed in a Functional-programming scheme (OOP not directly supported).
    * Currently no differential equations are provided, only the tools to solve them. The user must build their own diff-eqs in a separate script and then import the solvers from TidalPy.integration.
    * OOP classes can be used in the solvers but must be wrapped by njit-capable functions.
        * For example: Using a planet core's mass that was calculated via Burnman can be used like:
~~~~
from TidalPy.planets import build_burnman_world
from TidalPy.integration.multi_rheo_integration import multi_rheo_integration
from numba import njit

planet = build_burnman_world('charon')
core_mass = planet.core.mass

@njit
def my_custom_diffeq(time, dependent_variables):
    # Do stuff
    core_density = core_mass / core_volume
    # Do stuff with core_mass (referenced as a global variable)
    return change_in_dependent_variables	
~~~~

* Major Changes
    * Replaced the custom logging system with one based on Python's logging package.
    * Modified OOP Backend
        * Models are now free to access the state properties of layers and planets. This eliminates the need to pass arguments into model classes for calculations. All the user has to do is change the layer and/or planet's state properties and then those changes will automatically propagate.
    * New Rheology Scheme
        * To better follow real physics, all strength-based models have been moved under the redesigned `Rheology` class.
        * This includes: liquid and solid viscosity calculations, partial melting, and complex compliance
        * Love numbers are now calculated by a new `Tides` class instead of the `Rheology` (see next bullet point).
    * New Tides Module
        * New `Tides` class and child classes have been implemented which handle all Love number and tidal calculations for both a rheology-based approach or for the CPL/CTL model.
        * This module contains functionality to calculate tidal heating and tidal potential derivatives based on a complex compliance function and the various thermal and orbital parameters.
        * The functions to calculate the Love numbers are also stored here.
        * New eccentricity functions have been added including terms up to and including e^20 and tidal order l=7.
        * New inclination functions have been added (with arbitrary accuracy) up to tidal order l=7.
    * Other
        * Various functions and classes have been reworked or removed to simplify the code base. 

* Minor Changes
    * Removed ".pyname" from classes.
    * Utilities package has been reorganized.
    * Changed the setup pipeline to only require one command.
    
* QOL Improvements
    * Many new docstrings, type hints, and overall clean up of functions and classes. All docstrings now follow the numpy format.
    * Many new tests for both the functional and OOP versions of TidalPy.
    * More comments and spelling/typo fixes everywhere.
    
## Version 0.1.0 Alpha (July 2019)
Main Release - Changes not tracked