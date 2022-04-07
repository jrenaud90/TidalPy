# TidalPy Major Change Log

### Version 0.3.5 Alpha (Spring 2022)

* Minor Changes
  * Removed the soon-to-be deprecated np.float, np.complex, np.int references.
  * Added a check to see if cartopy is installed before functions that depend on it are imported.
* Bug Fixes
  * Fixed coverage problem with github actions.
  * Fixed issue importing nbTuple.
  * Fixed issues with GitHub actions and coverage.

### Version 0.3.4 Alpha (Winter/Spring 2022)

*Multilayer scripts based on 0.3.3 or earlier will likely break with this version!*

* Major Changes
    * Added `TidalPy.modes.multilayer_modes.py` module to offer simplified calculation of multilayer tidal
      heating.
    * Added `GridPlot` class to quickly make grid-like matplotlib figures.
      Checkout `TidalPy.utilities.graphics.grid_plot.py`
    * Added `Cartopy` dependence.
        * Can now make cool projection maps! Added basic functionality to `TidalPy.utilities.graphics.global_map.py`.
        * New jupyter notebooks to showcase map projects and GridPlot functionality.
    * Improved performance on both mode and non-mode tidal potential functions by at least a factor of 3. If used
      correctly these can be nearly 100x faster.
    * Added a new obliquity version of the mode version tidal potential.
    * Stress and strain relationship for multi-layer tides now allows for arbitrary rheology.
    * Created a single multilayer solver to handle an arbitrary layer structure.
      See `TidalPy.tides.multilayer.numerical_int.solver.py`
    * Stress & Strain relationship now accounts for arbitrary rheology.
    * Created a single multi-mode solver for multilayer problems. See `TidalPy.tides.modes.multilayer_modes.py`
    * TidalPy now defaults to using the frequency dependent zeta versions of Andrade and Sundberg rheologies.
        * This was done to avoid issues with real(complex_comp) at zero frequency which happens in multi mode
          calculations.
    * Added numba-safe version of multilayer calc

* Minor Changes
    * Added a helper function to quickly calculate masses, volumes, and gravity for spherical shells provided a radius
      and density array: `TidalPy.utilities.spherical_helper.calculate_mass_gravity_arrays`.
    * Added more colormaps, updated how reserved versions are constructed, updated old maps.
    * Better support for post-multiprocessing function inputs.
    * Refactored tidal mode calculation functions into a new TidalPy.tides.modes module.
    * Added a voxel calculator to TidalPy.utilities.spherical_helper. Also, added related tests & benchmarking tools.
    * Added a dictionary to track known color maps. It can be imported at `TidalPy.utilities.cmaps.KNOWN_CMAPS`
    * Made some improvements to the unique_path function in io_helper.py
    * Rearranged the tidal potential argument order.
    * Added some sanity checks on the various kinds on both mode and non-mode tidal potentials to compare with one
      another.
    * Updated stress, strain, and displacement calculations to account for new low-memory calculation method.
    * Greatly increased performance of multilayer stress, strain, and potential calculations.
    * Refactored much of the multilayer functions from TidalPy.toolbox to TidalPy.tides.multilayer.numerical_int and sub
      modules
    * Reworked numba-safe RK integrator. Does not reproduce scipy exactly for chaotic functions when numba is on.

* Bug Fixes
  * Fixed issue that was causing some multilayer boundary calculations to be 10x slower.
  * Fixed issue where 2D arrays would not work on njited `neg_array_for_log_plot` function.
  * Fixed issue where there could be negative love numbers with frequency is zero in multi-mode calculator

### Version 0.3.3 Alpha (Fall 2021)

* Major Changes
    * Updated to work with BurnMan v1.0 release.
    * Updated the NSR tidal potential for multilayer code to a multi-modal version
    * Improvements to the multilayer module
* Minor Changes
    * Added more tests
    * Removed the remaining `_array` functions

### Version 0.3.2 Alpha (Summer 2021)

* Major Changes
    * Changes to setup to use released version of BurnMan.

### Version 0.3.1 Alpha (Summer 2021)

* Major Changes
    * Added a simplistic (not modal version) NSR tidal potential for use in the multilayer code.
    * Added multiprocessor toolkit to `TidalPy.utilities.multiprocessing`
* Minor Changes
    * Did large scale code reformatting (just style, no refactoring) on nearly all files.
    * Cleaned up some doc strings.
* Bug Fixes

## Version 0.3.0 Alpha (Spring 2021)

*Scripts based on 0.2.x will likely break with this version!*

* Major Changes:
    * Added the first iteration of a multilayer tidal calculator module in `TidalPy.tides.multilayer` this module
      provides basic functionality to calculate tidal dissipation in a semi-homogeneous, shell-based approach. This is
      more accurate than the pure homogeneous model used throughout the rest of TidalPy. The downside with the current
      version is that it does not allow for NSR or high eccentricity / obliquity. A future version will attempt to add
      in a more robust Tidal Potential equation which will allow for additional physics.
    * Setup.py has been revamped as has the installation process. This is in prep to allow for TidalPy to become
      available on PyPI.
    * Did away with most of the `_array` functions. Found a way for njit to compile a function to handle either arrays
      or floats.
        * Left the `self._func_array` (in addition to `self._func`) in the `model.py` classes just in case we ever **
          do** need to define array functions in the future: all the infrastructure is still in place.
    * Added a numba-safe Explicit Runge-Kutta integrator. This is fully wrapped in njit'd functions.
        * On its own this can be 5--20 times faster than `scipy.solve_ivp`.
        * This also allows the integration function to be used from within another njit'd function(s).
        * This is still very experimental so please use caution and testing when using it.
* Minor Changes
    * Numerous bug fixes.
    * Removed the array versions of the dynamic functions.
        * `use_array` is still tracked in OOP and some quick calculation functions. These may all be not necessary now.
    * New cookbook to showcase the multilayer calculations.
    * Added surface area slices to base physical class.
    * Fixed some issues with how radius slices are tracked within layers and worlds.
        * `<world/layer>.radii[0]` is never the radius at the bottom of the object, it is always one dx up.
        * An interpolation had to be added to Burnman layers since Burnman radii starts at 0.
    * Improved various docstrings.
    * Refactored the `TidalPy.tools` to `TidalPy.toolbox`.
    * Refactored `Cookbooks` to `Demos`.
    * conversions.semi_a2orbital_motion and orbital_motion2semi_a now always return np.nan where they used to return
      complex numbers.

## Version 0.2.1 Alpha (Fall 2020)

*Will break studies based on previous versions*

Note: TidalPy version of "0.2.0" was never made publicly available. This is the next version after 0.1.0.

* Major Changes
    * Major refactoring all over.
    * Replaced the custom logging system with one based on Python's logging package.
        * Logger will not print to console if using a Jupyter notebook.
        * User can decide if the log is saved to disk or not.
    * Modified OOP Backend
        * Models are now free to access the state properties of layers and planets. This eliminates the need to pass
          arguments into model classes for calculations. All the user has to do is change the layer and/or planet's
          state properties and then those changes will automatically propagate.
    * New Rheology Scheme
        * To better follow real physics, all strength-based models have been moved under the redesigned `Rheology`
          class.
        * This includes: liquid and solid viscosity calculations, partial melting, and complex compliance
        * Love numbers are now calculated by a new `Tides` class instead of the `Rheology` (see next bullet point).
    * New Tides Module
        * New `Tides` class and child classes have been implemented which handle all Love number and tidal calculations
          for both a rheology-based approach or for the CPL/CTL model.
        * This module contains functionality to calculate tidal heating and tidal potential derivatives based on a
          complex compliance function and the various thermal and orbital parameters.
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
    * Many new docstrings, type hints, and overall clean up of functions and classes. All docstrings should now follow
      the numpy format.
    * Many new tests for both the functional and OOP versions of TidalPy.
    * setup.py no longer requires a separate command line call to install the Burnman package.
    * More comments and spelling/typo fixes everywhere.
    * Added CVD-friendly color maps made by Crameri (2018; http://doi.org/10.5281/zenodo.1243862) to the
      utilities.graphics module
    * More log.debug() calls all over. This should hopefully help bugfixes in the future (especially OOP bugs).

* Bug Fixes:
    * Fixed bug in world_builder.py : build_from_world where nested dicts were not being overwritten as expected.
    * Fixed bug in melt fraction checks for the float version of the Henning model.
    * Fixed bug in the Henning melting model where viscosity and shear was not being calculated correctly during the
      breakdown band (critical melt fraction + ~5%). This only affected a phase space that was rarely important for
      tidal calculations (planets would pass through it *very* quickly).

* In the works:
    * Integration module
        * The beginning of an integration module was added. This is still very early and not particularly useful. Users
          are still recommended to build their own integration tools that wrap TidalPy functions.

## Version 0.1.0 Alpha (July 2019)

*Initial Release - Changes not tracked*
