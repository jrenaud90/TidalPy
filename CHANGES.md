# TidalPy Major Change Log

### Version 0.7.0 (2025-xx-xx)

#### New
- Added `TidalPy.utilities.math.optimize` module and associated tests.
  - `optimize.bracket` mimics the root finding algorithm of BurnMan's `bm.utils.math.bracket`.

### Version 0.6.9 (2025-09-19)

#### Fixes
* Fixed issue where `TidalPy.RadialSolver.shooting` would pick the incorrect starting index. If the starting layer (set by the starting radius) was not the first layer it could cause a int overflow and lead to access violation crashes. 

#### Dependencies
* Updates some GitHub action dependencies.

#### Tests
* Added more tests to check `RadialSolver` starting radius conditions to try to catch bugs like this patch fixed!

### Version 0.6.8 (2025-08-19)

#### Changes
* `TidalPy.tides.modes.collapse_multilayer_modes` now checks for ill-formed `radius_arrays` and will raise an error if it detects a problem.
* Improved how `projection_map` displays colorbar numbers.

#### Fixes
* Fixed issue where `TidalPy.tides.modes.collapse_multilayer_modes` was not treating macro-layer boundaries correctly.

### Version 0.6.7 (2025-08-19)

#### Changes
* `TidalPy.tides.modes.collapse_multilayer_modes` now requires a tuple of rheology instances for both shear and bulk rheology. One for each macro layer. This allows different rheologies to be used for different layers.

#### Fixes
* Fixed some warnings that were showing up in the multilayer demo notebook.

### Version 0.6.6 (2025-08-15)

#### Changes
* Updated to work with CyRK v0.15.1
  * This change results in about a 18% performance improvement for RadialSolver depending on the problem.
* Converted RadialSolver and EOSSolver status message to C++ strings.
* Updated some documentation.

#### Fixes
* Fixed potential issue where RadialSolver could use a very small or even negative maximum step size.

### Version 0.6.5 (2025-08-12)

#### New
* Added test to check if structure arrays have been changed.
* Added debug flag to installation files to help with cython debugging. 

#### Fixes
* Fixed issue where TidalPy structures (layers, planets, etc.) would return editable arrays instead of copies of arrays. This could lead to subsequent functions (like planet paint) changing the arrays. This fixes GitHub Issue [#74](https://github.com/jrenaud90/TidalPy/issues/74).
* Fixed issue where shooting method was corrupting data when starting radius was too large (which happened for higher orders of l) This fixes GitHub Issue [#72](https://github.com/jrenaud90/TidalPy/issues/72).
* Fixed broken demo notebooks.

### Version 0.6.4 (2025-04-10)

#### Dependencies
* Removed max version limit for platformdirs package

### Version 0.6.3 (2025-04-09)

#### Changes
* Changes GitHub actions to avoid testing when not needed.
* Now supports Python 3.13

#### Fixes
* Removed various files that were being included in the manifest

#### Dependencies
* Updated to work with CyRK 0.13.3
* Updated to work with Numpy 2.x

### Version 0.6.2 (2025-03-28)

* Bumping version to integrate into conda-forge.
* Fixed some build processes that could cause both numpy 1.X and 2.X to be installed and inconsistent use (only affected MacOS)

### Version 0.6.1 (2025-01-11)

#### Benchmarks and Performance
* Added performance benchmark for `TidalPy.RadialSolver.helpers`.

#### Fixes
* Fixed memory leak that was occurring in the EOS Solver (therefore also RadialSolver) due to CyRK's CySolverSolution not being dereferenced (this is a problem with CyRK, but hack applied to TidalPy until a fix can be made to CyRK).

#### Performance
* Improved `TidalPy.RadialSolver.helpers.build_rs_input_from_data` performance by factor of 10x to 150x depending on layer structure.
* Improved `TidalPy.RadialSolver.helpers.build_rs_input_homogeneous_layers` performance by factor of 15% to 3.5x depending on layer structure.
* Improved `TidalPy.RadialSolver.radial_solver` performance by around 10% depending on layer structure.

## Version 0.6.0 (2025-01-07)
#### Removed
* Removed support for the older non-cythonized `TidalPy.radial_solver` module in favor of `TidalPy.RadialSolver`

#### RadialSolver Changes
* Moved RadialSolver's Boundary Condition finder to its own function in `TidalPy.RadialSolver.boundaries.surface_bc.pyx` to allow it to be used by both the shooting and propagation matrix techniques.
* Decoupled radial solver from shooting method.
  * Moved the shooting method (formerly just called `cf_radial_solver`) to a dedicated file to prep for a different dedicated file for the prop matrix solver. 
  * Now `TidalPy.RadialSolver.solver` only contains driver functions and output structures.
* Added Propagation Matrix technique to RadialSolver
  * This is simplified for now. Only planets with 1 solid, static, incompressible layer are allowed. 
    * Other assumptions can be approximated, e.g., liquid layers use a small shear modulus.
    * Multiple layers should also work if you have discontinuities in density, shear, etc. within your "one layer".
  * A cythonized solid fundamental matrix implementation can be found in `TidalPy.RadialSolver.PropMatrix.solid_matrix`.
* New radial slicing scheme required:
  * TidalPy now requires `radial_solver` input arrays to be defined in a precise manner:
    * `radius_array` must start at 0.
    * Each layer's upper and lower radius must be in the `radius_array`. That means if there is more than one layer there will be two identical radius values!
      * E.g., if a planet has a ICB at 1000km and a CMB at 3500km. Then `radius_array` must be setup with 2 values of 1000km and 2 values of 3500km. 
      * Other parameters should be defined on a "as layer" basis. So shear modulus at the 1st 1000km would be the shear of the inner core, at the 2nd 1000km it would be the shear modulus of the outer core. Likewise shear modulus at the 1st 3500km would be for the outer core and at the 2nd 3500km would be the shear modulus for the mantle. Same goes for density and bulk modulus.
* Added warning to check for instabilities (based on large number of steps taken; requires `warnings=True`).
* Changes to `radial_solver` arguments:
  * Many changes to the order as well as additions and removals of arguments to `radial_solver` highly suggest looking through the updated documentation.
  * `radial_solver` no longer requires `gravity_array`.
    * New with this update is a self-consistent equation of state solver (EOSS) that is called from `radial_solver`. This EOSS is used to find gravity(r).
  * Bulk modulus must now be provided as a complex-valued array
    * If a non-zero imaginary portion is provided (e.g., found via a rheology) then bulk dissipation is now tracked.
    * If imaginary portions are zero, then bulk dissipation is ignored as in TidalPy v0.5.x and earlier. (this is actually dependent on your bulk rheology; it could also cause infinities...)
  * Several arguments have had slight name refactors which will break calls that used keyword arguments. Review the RadialSolver documentation or "TidalPy/RadialSolver/solver.pyx" for the new argument names.
  * `upper_radius_bylayer_array` must now be provided as a numpy array (previously a tuple of floats was acceptable).
  * Added new optional argument `surface_pressure` (default=0.0) used with EOSS to find pressure convergence.
  * Added new optional argument `core_model` (default=0) used to set the lower boundary condition when the propagation matrix technique is used.
  * Added new optional argument `starting_radius` (default=0.0) to allow the user to set the initial radius for the radial solution.
    * Setting the solution radius higher in the planet can improve solution stability when using the shooting method. Particularly if looking at higher harmonic `degree_l`s. There is a trade off with accuracy so advise testing.
    * The starting radius must be >= 0.0, if 0.0 is provided (the default) then TidalPy will use the Martens technique to find a suitable starting radius (function of `degree_l`, planet radius, and the new optional argument `start_radius_tolerance` which defaults to 1.0e-5).
  * Removed `limit_solution_to_radius` argument.
  * Added new optional argument `use_prop_matrix` (default=False) to use the propagation matrix technique over the shooting method.
  * Equation of State Solver arguments:
    * `eos_method_bylayer` - EOSS method to use for each layer (currently only "interpolate" is supported).
    * `eos_integration_method` Runge-Kutta method to use for EOSS (default="RK45"). `eos_rtol` and `eos_atol` can also be provided to control integration error.
    * `eos_pressure_tol` (default=1.0e-3) and `eos_max_iters` (default=40) control the pressure convergence of the EOSS.
  * Added optional argument `perform_checks` (default=True) that performs many checks on the user input before running the solution (small performance penalty, but highly recommend leaving on until your inputs are tested).
  * Added optional argument `log_info` (default=False) that will log key physical and diagnostic information to TidalPy's log (which can be set to be consol print, log file, or both via TidalPy's configurations). 
    * Note there is a performance hit when using this, particularly if logging to file is enabled.

**New RadialSolver Helpers**
* To assist with the generation of valid inputs to `radial_solver`, TidalPy now provides two helper functions:
  * For planets with homogeneous layers: `from TidalPy.RadialSolver import build_rs_input_homogeneous_layers` takes in attributes for a planet made of layers with constant physical properties and then provides the arrays and other required `radial_solver` inputs that conform to the new `radius_array` scheme.
  * For planets with inhomogeneous layers: `from TidalPy.RadialSolver import build_rs_input_from_data` which parses data arrays (loaded from a file or built elsewhere like using a more robust EOS than TidalPy offers) and makes changes to ensure they will work with `radial_solver`.

**Expanded RadialSolverSolution Attributes and Methods**
* The output of `radial_solver`, an instance of the `RadialSolverSolution`, has been greatly expanded to provide much more data and functionality to the user. Full details can be found in the new RadialSolver documentation. Highlights include:
  * EOSS results like `<solution>.mass`, `<solution>.moi`, `<solution>.moi_factor`, `<solution>.central_pressure`, `<solution>.surface_gravity`.
  * Diagnostic data like number of integration steps required per layer to find a solutions `<solution>.steps_taken`, or EOSS diagnostics: `<solution>.eos_iterations`, `<solution>.eos_pressure_error`, `<solution>.eos_success`, `<solution>.eos_message`, `<solution>.eos_steps_taken`.
    * Much of the new diagnostic data as well as key results can now be quickly printed using `<solution>.print_diagnostics()`.
  * Method to quickly plot the viscoelastic-gravitational solution y's `<solution>.plot_ys()`.
  * Method to quickly plot the EOS results `<solution>.plot_interior()`.
  * In addition to the previously provided attributes like `<solution>.love`, `<solution>.k`, `<solution>.h`, `<solution>.l`, `<solution>.result`.

#### Cython / C Changes
* Shifted away from `PyMem_Free` to `CyRK.utils.mem_free` to allow for consistency in future development.
  * Avoiding using these manual heap allocations whenever possible. Many new uses of smart pointers and C++ vectors.

#### Other Changes
* Updated GitHub actions.
* Moved to more consistent and robust types (e.g., int for degree_l vs. prior unsigned-char; unsigned-int).
* Added inverse function `cinv` in `TidalPy.utilities.math.complex`.
* New "TidalPy/constants.pyx" isolates all TidalPy constants. Available to both Python and Cython. Refactored all files to use the constants in this file.
* New numerics module `TidalPy.math.numerics` for low-level floating point functions.
  * New cythonized `isclose` function that matches functionality of python's `math.isclose`
* Cythonized radial sensitivity to shear/bulk functions in `TidalPy.tides.multilayer.sensitivity` (based on Tobie+2005)
* Cythonized radial heating functions that use the sensitivity to shear/bulk functions in `TidalPy.tides.multilayer.heating` (based on Tobie+2005)
* Improved logging so it is less spammy.
* Logger now logs all exceptions raised.
* Moved TidalPy's default config and world config dir to user's "Documents" folder (from system appdata folder). 
  * If upgrading from previous version of TidalPy, you can safely delete the old config directory.
    * On Windows the old dir was: "'C:\\Users\\<username>\\AppData\\Local\\TidalPy'"; The new dir is "'C:\\Users\\<username>\\Documents\\TidalPy'"
    * On Mac the old dir was: "'/Users/<username>/Library/Application Support/TidalPy'"; The new dir is "'/Users/<username>/Documents/TidalPy'"
    * On Linux the old dir was: "'/Users/<username>/.local/share/TidalPy'"; The new dir is "'/home/<username>/Documents/TidalPy'"
* New switch `TidalPy.log_to_file()` to quickly turn on saving log to file (this can also be adjusted in the TidalPy configurations).
* TidalPy now looks for an environment variable "TIDALPY_TEST_MODE" to turn on test mode during first initialization (can later be changed using the `TidalPy.test_mode()` command or setting `TidalPy._test_mode = False; TidalPy.reinit()`).
* Made use of more TidalPy-specific exceptions.
* Tweaked the `TidalPy.utilities.graphics -> yplot`.
* User can now override TidalPy.config using `TidalPy.reinit(<new config toml file path; or dictionary of configs>)`.
* Refactored and made improvements to `TidalPy.utilities.graphics.planet_plot`.

#### Dependencies
* Added support for CyRK v0.12.x
* Added support for Burnman v0.2.x

#### Documentation
* Reworked TidalPy's documentation structure in prep for a shift to Sphinx in the future.
* Greatly expanded and improved RadialSolver documentation which can be found in "TidalPy/Documentation/RadialSolver"

#### Fixes
* Fixed issue where `radial_solver` arrays could dealloc while references still pointed to them (hanging pointers).
* Missing Cython compile arguments in TidalPy's utilities, `nondimensional.pyx`.
* Fixed incorrect type in dynamic liquid layers that may have been causing some errors to propagate.
* Fixed issue where `TidalPy._config_path` was not being updated.
* Fixed issue where log files could not be created.

### Version 0.5.5 (2024-11-11)

Fixes:
* Fixed dependency compatibility issues.
* Fixed incorrect function signature type for scipy's `spherical_jn`. SciPy v.1.14.X uses a new signature which is breaking on MacOS. Limiting to "SciPy<1.14" for now. See GitHub Issue #65

### Version 0.5.4 (2024-04-30)

Fixes:
* Fixed RadialSolver frequency warning message for dynamic liquid layers not displaying for correct layer types.

Additions:
* Added way to suppress warning messages in RadialSolver.
  * To turn this suppression on, pass `warnings=False` to `RadialSolver.radial_solver`.

### Version 0.5.3 (2024-04-26)

Fixes:
* RadialSolver: Fixed bug where solutions between liquid and solid layers were not propagating correctly.

Additions:
* New Love number benchmarks for Earth provided by [Nick Wagner](https://github.com/nlwagner) in `Benchmarks & Performance\RadialSolver\Earth Love Numbers.ipynb` (Jupyter Notebook).

Changes:
* Pre-allocated several cythonized arrays to nans to help with debugging.
* Provided more error messages to improve user experience.
* Cythonized non-dim function now takes in the planet's density and radius as variables to change.
* Improved the Tobie and Roberts benchmarks for radial solver.

Other:
* Updated to work with CyRK 0.8.7

### Version 0.5.2 (2024-02-22)

Documentation
* Improved RadialSolver documentation regarding higher degree-l.
* Added info about issues that can arise from using non C-continguous arrays in cythonized functions.

Fixes:
* Added error message to `RadialSolver.radial_solver` if length of provided assumption tuples is not the same.
* Fixed issue where non C-continguous arrays were allowed in cythonized functions that required them.

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
