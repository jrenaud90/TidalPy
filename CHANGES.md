# TidalPy Major Change Log

### Version 0.6.0 *PLANNED*
* Other Major Changes
  * Remove support for the older non-cythonized `radial_solver` module.

### Version 0.5.1 (2024-02-14)
* Removed Python 3.8 support due to issues with building SciPy.

### Version 0.5.0 (2024-02-14)
_This version is likely to break code based on TidalPy v0.4.X and earlier_

Cythonizing TidalPy
* A major change starting with v0.5.0 is the switch from numba.njited functions to cython precompiled functions and
extension classes. The reasons for doing this are numerous. This transition will be completed in stages
with minor versions (v0.X.0) each bringing a new set of cythonized updates until all njited functions are retired.
* For this version: 
  * Converted `TidalPy.radial_solver.radial_solver` to cythonized `TidalPy.RadialSolver.radial_solver`.
    * The old radial solver method will be removed in TidalPy version 0.6.0.
  * Added new cython-based `TidalPy.utilities.classes.base_x` base cython extension class that other classes are built off of.
  * Converted `TidalPy.rheology.complex_compliances.compliance_models` to cythonized `TidalPy.rheology.models`.
    * Improved the new rheology methods to better handle extreme values of frequency and modulus.
    * The old rheology solvers will be removed in a future release of python.
  * Added several new cython-based helper functions in the utilities package.

Other Major Changes
* Added support for Python 3.11 and 3.12. TidalPy now runs on Python 3.8--3.12.
  * Note that currently the Burnman package does not support 3.12 so that functionality is limited to python 3.8-3.11.
* Removed support for `solver_numba` in the `radial_solver` module.
* Removed some imports from main package and sub modules' `__init__` to avoid slow load times.
* Moved conversion tools from `TidalPy.toolbox.conversions` to `TidalPy.utilities.conversions`.
* Changed setup files so that cython code can be compiled.
  * `special` - for high-performance, general, scientific functions.
* Moved TidalPy configs to a standard user directory by default. The specific location will depend on the OS.
  * Default configs will be installed on the first `import TidalPy` call after installation.
    * These defaults are stored in the `TidalPy.defaultc.py` as a string which is copy and pasted to the new `TidalPy_Configs.toml`.
  * There is a new `TidalPy.clear_data()` function to delete all data stored in these locations. Data will be rebuilt the next time TidalPy is imported.
  * New `TidalPy.set_config(config_path)` to change the active configuration file used by TidalPy.
    * Note that `TidalPy.reinit()` should be called after changing the configurations.
  * New `TidalPy.set_world_dir(world_dir_path)` to change which directory to pull world configs from. 
  * Moved away from the system of `default.py` configurations for sub modules. All default configs are stored in the same `TidalPy_Config.toml`
* Shifted from `json` to `toml` files for world configs.
  * Store all world configs to a zip file for easier distribution.
* TidalPy now requires:
  * CyRK>=0.8.6
  * Cython>=3.0.0
* Moved `BurnMan` 3rd party dependence to a more dedicated `Extending` folder for future development.
* To make TidalPy lighter weight we are starting to remove a lot of 3rd party packages.

Minor Changes and New Features
* `complex_compliance` configurations are now stored in the top level `rheology` in all configs.
  * For example, in prior TidalPy versions you would need to change the complex compliance model by editing `config['layers']['mantle']['rheology']['complex_compliance']['model'] = 'andrade'`. Now this would be: `config['layers']['mantle']['rheology']['model'] = 'andrade'`.
* Added unique frequency finder functions to the `modes` module.
* Moved most of the type hints behind the `typing.TYPE_CHECKING` flag.
* Moved non-critical files out of repository.
* Created a new `tides.heating` module and moved the volumetric heating calculations there.
* Expanded the performance suite to better track the `radial_solver` module.
* Moved `cache.py` to top-level.
* Turned off numba cacheing on several functions that may be used in the radial solver.
  * rheology
    * complex compliance functions
  * radial_solver.numerical
    * initial guess functions
    * interface functions
* Converted radial_solver.numerical initial guess and interface functions output to np.ndarrays rather than numba lists.
* Removed `config_helper.py` and the functions defined within.
* New RadialSolver class now supports more than just boolean inputs.
  * Future proofing to allow for a greater variety of layer types.
* Added exoplanet archive download functionality in `TidalPy.utilities.exoplanets`.
  
Bug Fixes
* Fixed floating point comparison bug in `multilayer_modes` solver.
* Fixed obliquity not being used issue in quick tides calculator.
* Fixed issue in incorrect TidalPy version being loaded into the package.

Performance Improvements
* Improved the performance of the stress and strain calculator by ~20%.
* Cythonize Performance Increases:
  * New `RadialSolver.radial_solver` leads to a ~50x performance boost.
  * New cythonized rheology models are 500% faster for arrays; 40,000% faster for scalars (not a typo!)

#### Version 0.4.1 Alpha (Spring 2023)
Major Changes
* Moved `radial_solver` to a top-level module of TidalPy.
  * Added `find_love` function to the `radial_solver` module for the calculation of Love and Shida numbers.
  * Added `sensitivity_to_shear` function to the `radial_solver` module based on Tobie et al. (2005).
  * Added `sensitivity_to_bulk` function to the `radial_solver` module based on Tobie et al. (2005).
* Added `newton` and `elastic` complex compliances to rheology module for ALMA comparisons.

Minor Changes
* Updated and added to `radial_solver` test suite.
* Changed `radial_solver.interfaces` functions to only require the top most value of the radial solutions rather than the whole array.

Bug Fixes
* Fixed bug in `radial_solver` where interfaces between a dynamic and static liquid were not correctly handled.

### Version 0.4.0 Alpha (Winter 2022/2023)

Major Changes 
* Added a multilayer tidal potential that allows for arbitrary obliquity.
* Added in load Love number calculations to the multilayer code.
* Removed a lot of 3rd-party dependencies to make TidalPy's install more lean.
* Switched over to using the integrators from the new `CyRK` package
  * Changed the signature of the numerical-int multilayer solver.
    * **Breaks Old Code**
* Issue with Numba 0.55 and dictionary updates. This restricts TidalPy to Python version 3.9 or lower.
* Started a lot of prep work for a move to sphinx or similar for documentation (this will be a 0.5.0 feature).
* Created a generalized multi-layer collapse function in `TidalPy.tides.multilayer.numerical_int.collapse`.
  * This will likely break old code. It was introduced in v0.4.0.dev10.
  * Arbitrary layer structures can now be fed into the y-solver. e.g., liquid-solid, solid-liquid, liquid-solid-liquid-solid, etc.
    * Note: Multiple dynamic liquid layers will likely lead to stability problems.
  * y-solver no longer requires a `model_name` variable. The function will automatically utilize the correct model based on the `is_layer_solid` and `is_layer_static` lists.
  * Removed the old code that handled individual cases.
* Interfaces between multiple liquid layers are now supported (but not well tested).
  * Two liquid layers of different types can be next to one another (dynamic-static; static-dynamic).
* y-solvers functional argument signature has changed (both regular and numba versions)
  * **Breaks Old Code**
  * Numerical integrator declaration has changed.
  * Planet_bulk_density is now a required argument.
  * Modified several other TidalPy functions to match this new call signature for the y-solver.
* Removed the y-derivative return from the propagation matrix (to match the output format of the numerical int)
  * **Breaks Old Code**
* Refactored `tidal_y_solver` to `radial_solver` since non-tidal calculations can be made with it.
  * **Breaks Old Code (pre v0.4.0.dev11)**
* Switched from `setup.py` to a streamlined `pyproject.toml` installation process.
* Changes to radial ODE's
  * Input arguments and output diffeqs are now passed as numpy arrays rather than tuples.
  * Input and outputs are now passed as floats not complex (doubling the number of terms)
* Added `numba-scipy` dependence to allow the use of scipy's special functions. 
  * Removed the pre-calculated factorial method. Using scipy's gamma now.
  * TODO: Note the numba-scipy package on github is not updated to the newest version of scipy. Packaging numba-scipy with TidalPy for now.

Performance Improvements
* Improved the performance of the pure-numba radial solve by ~10%

Minor Changes
* Added newer functions to the performance recording suite.
* Improved performance of eccentricity and inclination functions for non-multilayer tidal calculations.
* Added TidalPy to Zenodo. DOI added to readme.
* Removed Gitter account for now.
* Added files and functions to quickly install additional 3rd party applications
* Improvements to GitHub workflows
* Switched over to using the 3rd party `cmcrameri` colormap package rather than trying to maintain it within tidalpy
* Added a delta time to HH:MM:SS converter to utilities.string_helper
* Made improvements to multiprocessing user info
* Cleaned up & added some docstrings and type hints
* Tidal y solver for the prop matrix method no longer returns the y-derivatives.
  * dy1/dr is now calculated directly in the `decompose()` function.
* Created a config helper function `TidalPy.test_mode()` to quickly setup TidalPy configs for pytest'ing
* Cleaned up comments and reordered items in `multilayer.numerical_int.collapse`.
* Fixed a bug when using SciPy's Radau integrator method.
* The version number is now checked with importlib in TidalPy.__init__. Version number should only be changed in the pyproject.toml.
* Updated the pure-numba version of the radial solver's argument signature to better match the python implementation.
* Updated Github Actions

Bug Fixes
* Fixed bug in GridPlot related to number of subplots.
* Fixed bug in global variable for world config loader.
* Fixed type hint bug in the numba-based tidal y solver.
* Fixed bug that caused numba-based tidal y solver to not compile.
* Fixed bug that was causing full TidalPy log to print while in a Jupyter Notebook environment.
  * If you would like the log to print in a notebook then use `TidalPy.toggle_log_print_in_jupyter()` or set the
`print_log_in_jupyter` to `True` in the "configurations.py" file.
* Fixed an error in the multimode volumetric heating calculation where "_rr", "_thth", "_phiphi" were being double
counted
* Fixed an error in the stress/strain calculations where the static, instead of complex, shear was being used.

### Version 0.3.5 Alpha (Spring 2022)

Minor Changes
* Removed the soon-to-be deprecated np.float, np.complex, np.int references.
* Added a check to see if cartopy is installed before functions that depend on it are imported. 

Bug Fixes
* Fixed coverage problem with github actions.
* Fixed issue importing nbTuple.
* Fixed issues with GitHub actions and coverage.

### Version 0.3.4 Alpha (Winter/Spring 2022)

*Multilayer scripts based on 0.3.3 or earlier will likely break with this version!*

Major Changes
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

Minor Changes
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

Bug Fixes
* Fixed issue that was causing some multilayer boundary calculations to be 10x slower.
* Fixed issue where 2D arrays would not work on njited `neg_array_for_log_plot` function.
* Fixed issue where there could be negative love numbers with frequency is zero in multi-mode calculator

### Version 0.3.3 Alpha (Fall 2021)

Major Changes
* Updated to work with BurnMan v1.0 release.
* Updated the NSR tidal potential for multilayer code to a multi-modal version
* Improvements to the multilayer module

Minor Changes
* Added more tests
* Removed the remaining `_array` functions

### Version 0.3.2 Alpha (Summer 2021)

Major Changes
* Changes to setup to use released version of BurnMan.

### Version 0.3.1 Alpha (Summer 2021)

Major Changes
* Added a simplistic (not modal version) NSR tidal potential for use in the multilayer code.
* Added multiprocessor toolkit to `TidalPy.utilities.multiprocessing`

Minor Changes
* Did large scale code reformatting (just style, no refactoring) on nearly all files.
* Cleaned up some doc strings.

## Version 0.3.0 Alpha (Spring 2021)

*Scripts based on 0.2.x will likely break with this version!*

Major Changes:
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

Minor Changes
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

Major Changes
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

Minor Changes
* Removed ".pyname" from classes.
* Utilities package has been reorganized.
* Changed the setup pipeline to only require one command.

QOL Improvements
* Many new docstrings, type hints, and overall clean up of functions and classes. All docstrings should now follow
  the numpy format.
* Many new tests for both the functional and OOP versions of TidalPy.
* setup.py no longer requires a separate command line call to install the Burnman package.
* More comments and spelling/typo fixes everywhere.
* Added CVD-friendly color maps made by Crameri (2018; http://doi.org/10.5281/zenodo.1243862) to the
  utilities.graphics module
* More log.debug() calls all over. This should hopefully help bugfixes in the future (especially OOP bugs).

Bug Fixes:
* Fixed bug in world_builder.py : build_from_world where nested dicts were not being overwritten as expected.
* Fixed bug in melt fraction checks for the float version of the Henning model.
* Fixed bug in the Henning melting model where viscosity and shear was not being calculated correctly during the
  breakdown band (critical melt fraction + ~5%). This only affected a phase space that was rarely important for
  tidal calculations (planets would pass through it *very* quickly).

## Version 0.1.0 Alpha (July 2019)

*Initial Release - Changes not tracked*
