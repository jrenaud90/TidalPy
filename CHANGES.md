# TidalPy Major Change Log

## Version 0.2.1 Alpha (Fall 2020)
*Will break studies based on previous versions*

Note: TidalPy version of "0.2.0" was never made publicly available. This is the next version after 0.1.0.

* Major Changes
    * Major refactoring all over.
    * Replaced the custom logging system with one based on Python's logging package.
        * Logger will not print to console if using a Jupyter notebook.
        * User can decide if the log is saved to disk or not.
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
    * Many new docstrings, type hints, and overall clean up of functions and classes. All docstrings should now follow the numpy format.
    * Many new tests for both the functional and OOP versions of TidalPy.
    * setup.py no longer requires a separate command line call to install the Burnman package.
    * More comments and spelling/typo fixes everywhere.
    * Added CVD-friendly color maps made by Crameri (2018; http://doi.org/10.5281/zenodo.1243862) to the utilities.graphics module
    * More log.debug() calls all over. This should hopefully help bugfixes in the future (especially OOP bugs).
    
* Bug Fixes:
    * Fixed bug in world_builder.py : build_from_world where nested dicts were not being overwritten as expected.
    * Fixed bug in melt fraction checks for the float version of the Henning model.
    * Fixed bug in the Henning melting model where viscosity and shear was not being calculated correctly during the breakdown band (critical melt fraction + ~5%). This only affected a phase space that was rarely important for tidal calculations (planets would pass through it *very* quickly).

* In the works:
    * Integration module
        * The beginning of an integration module was added. This is still very early and not particularly useful. Users are still recommended to build their own integration tools that wrap TidalPy functions.

## Version 0.1.0 Alpha (July 2019)
Main Release - Changes not tracked