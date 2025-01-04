# Calculating Love Numbers using TidalPy's RadialSolver module and functions
_In addition to the documentation in this file, there is a demo notebook "3 - Calculate Love Numbers.ipynb" in the Demos folder. As well as notebooks in "Benchmarks / Radial Solver" folder that may be helpful._

TidalPy's `RadialSolver` package allows a user to estimate a planet's global, viscoelastic
[Love numbers](https://en.wikipedia.org/wiki/Love_number). These numbers can then be used to determine the magnitude of tidal
dissipation, speed of rotational/orbital changes, and provide predictions for gravity or distance measurements. 

TidalPy allows for both the use of a numerical shooting method and a propagation matrix technique. However, the 
shooting method is much further developed within TidalPy, has more options and control, and is the default method.
This method integrates a set of 3x6 viscoelastic-gravitational ordinary differential equations
from the core of the planet to its surface. The solution is calculated up to layer interfaces (layers are defined by
phase changes within the planet) and then new ODEs are solved within the next layer, repeating up to the planet's
surface. The final result is a super position of the different solutions in each layer, where each
solutions' coefficients are determined by boundary conditions at the surface of the planet. 

## References
To learn more about TidalPy's underlying methods please review these references.

**Shooting Method**:
- Takeuchi, H., & Saito, M. (1972). Seismic Surface Waves. In Methods in Computational Physics: Advances in Research and Applications (Vol. 11, pp. 217–295). Elsevier. https://doi.org/10.1016/B978-0-12-460811-5.50010-6
- Tobie, G., Mocquet, A., & Sotin, C. (2005). Tidal dissipation within large icy satellites: Applications to Europa and Titan. Icarus, 177(2), 534–549. https://doi.org/10.1016/j.icarus.2005.04.006
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., & Běhounková, M. (2021). Solid tides in Io’s partially molten interior: Contribution of bulk dissipation. Astronomy & Astrophysics, 650, A72. https://doi.org/10.1051/0004-6361/202039433

**Interfaces, Constants, and Assumptions**:
- Saito, M. (1974). Some problems of static deformation of the earth. Journal of Physics of the Earth, 22(1), 123–140. https://doi.org/10.4294/jpe1952.22.123
- Beuthe, M. (2015). Tidal Love numbers of membrane worlds: Europa, Titan, and Co. Icarus, 258, 239–266. https://doi.org/10.1016/j.icarus.2015.06.008

**Starting Conditions**:
- Kamata, S., Matsuyama, I., & Nimmo, F. (2015). Tidal resonance in icy satellites with subsurface oceans. Journal of Geophysical Research (Planets), 120, 1528–1542. https://doi.org/10.1002/2015JE004821
- Martens, H. R. (2016). USING EARTH DEFORMATION CAUSED BY SURFACE MASS LOADING TO CONSTRAIN THE ELASTIC STRUCTURE OF THE CRUST AND MANTLE. CalTech, PHD Thesis.

**Propagation Matrix Method**:
- Sabadini, R., & Vermeersen, B. (2004). Global dynamics of the earth: Applications of normal mode relaxation theory to solid-earth geophysics. Kluwer Academic Publishers.
- Roberts, J. H., & Nimmo, F. (2008). Tidal heating and the long-term stability of a subsurface ocean on Enceladus. Icarus, 194(2), 675–689. https://doi.org/10.1016/j.icarus.2007.11.010
- Henning, W. G., & Hurford, T. (2014). TIDAL HEATING IN MULTILAYERED TERRESTRIAL EXOPLANETS. The Astrophysical Journal, 789(1), 30. https://doi.org/10.1088/0004-637X/789/1/30

## Radial Solver Function `TidalPy.RadialSolver.radial_solver`
The `radial_solver` function, contained in the `TidalPy.RadialSolver` module is the main way to solve the radial 
functions from Python. There are also cython hooks for faster performance (see "4 - Cython API" documentation).

Notes:
- All arrays must be [C-contiguous](https://stackoverflow.com/questions/26998223/what-is-the-difference-between-contiguous-and-non-contiguous-arrays). If you suspect that an array may not be C-contiguous you can use the numpy function `arr = np.ascontiguousarray(arr)` to ensure that they are before being passed to the rheology methods.
- At least 5 slices per layer is required (so total size of arrays must be at least 5x num_layers).
- RadialSolver will solve an equation of state to determine various other required properties (such as gravity)
    - Currently, only an interpolation EOS is implemented, meaning that if properties change with radius (e.g., density(r)) within layers then the arrays must be robust enough to capture those changes.
    - For example, for a planet with homogenous layers then 5 slices per layer is plenty because density does not change within the layer.
    - On the other hand, if density(r) != constant within the layer then you need to ensure there are enough slices to capture these changes. This is even more important for parameters that change faster, like viscosity.
- The radius array (and all properties that change with radius) must follow this format:
    - Starts at r=0.0
    - Has r values at the top and bottom of each interface. Meaning if there are 2+ layers then the interface radius value will be in the radius array _twice_.
- TidalPy provides helper functions to easily create or modify arrays that follow RadialSolver's requirements. Please see "3 - Radial Solver Helpers.md" documentation for more details.

```python
from TidalPy.RadialSolver import radial_solver
rs_solution = radial_solver(
    # # # Required inputs # # #

    radius_array,
    # [m] (type: double array; size = total_slices)
    
    density_array,
    # [kg m-3] (type: double array; size = total_slices)
    
    complex_bulk_modulus_array,
    # [Pa] (type: double complex array; size = total_slices)
   
    complex_shear_modulus_array, 
    # [Pa] (double complex array; size = total_slices)

    # Note that `complex_shear_modulus_array` and `complex_bulk_modulus_array` are complex-valued array.
    # This is a result of applying It is the result of applying a rheological function to the respective modulus and viscosity.
    # The complex part of either array can be set to 0 if you do not care about shear / bulk dissipation.
    
    frequency,
    # Scalar forcing frequency [rad s-1] (type: double)
    
    planet_bulk_density,
    # Scalar bulk density of planet [kg m-3] (type: double)
    # The following tuples define the assumptions used for the major layers within a planet.
    
    layer_types,  
    # Tuple of layer types (type: tuple of strings)
    # Options:
    #  - 'solid'
    #  - 'liquid'
    # Propagation matrix only allows for a single solid layer.
    
    is_static_bylayer,
    # Is each layer using the quasi-static tidal assumption (type: tuple of bools)
    # See Beuthe (2015) for a detailed description of this assumption.
    # Propagation matrix only allows for a static layer.
    
    is_incompressible_bylayer,
    # Is each layer using the incompressible tidal assumption (type: tuple of bools)
    # Propagation matrix only allows for a incompressible layer.
    
    upper_radius_bylayer_array,
    # Upper radius of each layer [m] (type: double array; size = num_layers)
    

    # # # Below are optional arguments. The default values are shown after the "=". # # #

    degree_l = 2,
    # Harmonic degree used. 
    # Note that the stability of the solution gets worse with higher-l. You may need to increase the relative tolerance
    # or use simpler layer assumptions to get a successful solution. 
    # It particularly starts to break down for l > 10 but stable solutions exist up to l=40.
    # For higher degrees, it is recommended to start integration higher in the planet. I.,e. above the
    #  core-mantle-boundary (if applicable).
    
    solve_for = None,
    # What to solve for (type: tuple[str, ...])
    # Options:
    #  - 'tidal'    Tidal Love numbers
    #  - 'loading'  Loading Love numbers
    #  - 'free'     Numbers for a free surface condition
    # If set to "None" then "('tidal',)" will be assumed.
    # You must provide these as a tuple of strings, even if you are only solving for one thing at a time. 
    # For example: solve_for = ('tidal',)  <-- Note the required comma.
    # It is more computationally efficient to solve for multiple Love number types at same time if you need more.
    # For example: solve_for = ('tidal', 'loading')
    # The number of items passed to this variable sets the size of `num_solve_for` discussed later in this document.

    starting_radius = 0.0,
    # What radius should the solution start at? [m] (type: scalar double)
    # Some problems are more stable if the starting radius is higher in the planet. This tends to be the case 
    # when degree_l >> 2. If set to 0.0, the default, then TidalPy will use Martens (2016) technique to determine 
    # a good starting radius depending on the degree_l and the tolerance set in the next variable.

    start_radius_tolerance = 1.0e-5
    # Starting radius tolerance (type: scalar double)
    # Used to automatically determine a decent starting radius, see previous argument for more details.

    nondimensionalize = True,
    # If True, then the inputs will be non-dimensionalized before integration. Radial solutions will be
    # re-dimensionalized before solution is given back to the user.
    # Generally setting this to True will provide more solution stability.
    

    # # # The following indented options are only used for the shooting method # # #

    use_kamata = False,
    # Use starting conditions as defined in Kamata et al. (2015; JGR:P)
    # If set to False then the starting conditions from Takeuchi & Saito (1972) will be used
    #
    # Kamata's starting conditions are more stable for incompressible layers. In fact, for an incompressible solid core
    # at the center of the planet you must set `use_kamata=True` as this starting condition is undefined in TS72.
    
    integration_method = 'DOP853',
    # Integration method to use during integration
    # Options:
    #  - 'RK23'    Explicit Runge-Kutta method of order 3(2)
    #  - 'RK45'    Explicit Runge-Kutta method of order 5(4)
    #  - 'DOP853'  Explicit Runge-Kutta method of order 8
    
    integration_rtol = 1.0e-5,
    # Integration relative tolerance

    integration_atol = 1.0e-8,
    # Integration absolute tolerance (tolerance when y is near 0)

    scale_rtols_by_layer_type = False,
    # If True, then relative tolerance will be modified by the layer type (generally liquid layers require smaller rtol)
    # !! This is an experimental feature and has not been thoroughly tested.

    max_num_steps = 500_000,
    # Maximum number of steps allowed for each integration (note that multiple integrations can occur depending on the
    # number and type of layers). For most cases, the integration will only require <5000 steps per solution.
    
    expected_size = 1000,
    # Expected number of steps needed to complete adaptive radial integration. It is better to overshoot this.
    
    max_ram_MB = 500,
    # Maximum amount of ram the integrator is allowed. Note that real ram usage will be larger than this.
    
    max_step = 0,
    # Maximum step size allowed. If 0 then the integrator will attempt to find a reasonable maximum step size.
    
    # # # The following arguments are only used by the propagation matrix method # # #

    use_prop_matrix = False,
    # Tells RadialSolver to use the propagation matrix technique instead of the shooting method.
    # 
    # Note that this method is less accurate and has many more restrictions. It is included here for comparisons purposes.
    # 
    # Note, unlike the shooting method, the matrix method is much more sensitive to the size of the input arrays (`total_size`).
    # This size sets the dimension of the propagation matrix which, if too small, will lead to large errors.
    
    core_model = 0,
    # Sets the inner core assumptions used by the propagation matrix technique. Ignored by the shooting method.
    # 
    # If the starting radius is near the planet's core then the final solution is not very sensitive to the starting
    # condition. However, if you are starting higher up in the planet it becomes more important and you may want to 
    # test multiple versions.
    # 
    # Options:
    #     * 0: Henning & Hurford (2014), uses seed matrix from propagation matrix
    #     * 1: Roberts & Nimmo (2008), very small liquid core (basically a solid core but with a tiny liquid core)
    #     * 2: Henning & Hurford (2014), Solid inner core
    #     * 3: Tobie+ 2005, liquid inner core (as determined by Marc Neveu for IcyDwarf)
    #     * 4: Sabadini & Veermeerson (2004), More complex interface matrix


    # # # The following arguments are for TidalPy's equation of state solver # # #

    eos_method_bylayer = None,
    # Tuple of EOS methods for each layer. This is a tuple of strings.
    #    If `None` then will use the default for each layer (interpolation)
    #    Currently supported methods:
    #        - "interpolation"

    surface_pressure = 0.0,
    # Pressure at planet's surface [Pa] (type: scalar double)
    # Used by the equation of state solver to find the interior pressure throughout the planet.

    eos_integration_method = 'DOP853',
    # Integration method used to solve for the planet's equation of state (type: str)
    # Same options as `integration_method` discussed above.

    eos_rtol = 1.0e-3,
    # Integration relative tolerance for equation of state solver.

    eos_atol = 1.0e-5,
    # Integration absolute tolerance for equation of state solver.

    eos_pressure_tol = 1.0e-3,
    # Tolerance used when fitting to the surface pressure in the equation of state solver. (type: scalar double)

    eos_max_iters = 40,
    # Maximum number of iterations used to converge surface pressure in equation of state solver.


    # # # The following arguments are for record keeping and debugging # # #

    verbose = False,
    # If True, then additional information will be printed to the terminal during the solution. 

    warnings = True,
    # If True, then warnings will be printed to the terminal during the solution if they arise.
    
    raise_on_fail = False,
    # If True, then the solver will raise an exception if integration was not successful. By default RadialSolver fails silently.

    perform_checks = True,
    # Performs sanity checks that raise python exceptions. If turned off then these checks will be skipped providing 
    # some boost to performance but at the risk of uncaught exceptions (crashes).

    log_info = False
    # Flag to turn on logging of key information (diagnostic and physical) from the RadialSolverSolution.
    # Note there is a performance hit if this is true, particularly if file logging is enabled.
    )
```

## `RadialSolverSolution` Class
The output of `radial_solver` is a `RadialSolverSolution` class. This class has several helpful attributes described below.

In addition to the full `y` radial solution results, you can also quickly access the Love numbers using:


**Performance Note**

There is a minor overhead when accessing any of the Solver's attributes (e.g., `.result; .love; .k; .h; .l`).
If you have a code that accesses these numbers often (more than once) it is better to store them in a local variable.
However, please note the common problem discussed in the "Undefined or random values for Love numbers" section.

```python
k_local = solution.k  # Performs a background lookup and np.ndarray operation to produce an array for all k Love numbers.

# ... Do a bunch of stuff with the new "k_local" variable.
```

### Radial Solution Troubleshooting

#### Undefined or random values for Love numbers
This issue is further described and discussed in GitHub Issue [#59](https://github.com/jrenaud90/TidalPy/issues/59). 
A common problem when using radial solutions is losing the variable's scope and the values of, for example, the Love
numbers becoming undefined. 

For example, consider this wrapper function:
```python
def wrap():
    result = radial_solver(...)
    k_inner = result.k # k_inner is a numpy array of double complex floating point values.
    return k_inner

k_outer = wrap()
```
TidalPy currently stores results in a heap allocated array which is returned when `result.k` is invoked. However, when
`result` (which is a `RadialSolverSolution` class instance) loses reference at the end of the `wrap` function then this
array is deallocated. The result is that `k_outer` now points to unallocated memory instead of the desired values.

To work around this, make a copy of the array whose memory is then handled by the python interpreter and not deallocated
when the scope changes:

```python
import numpy as np

def wrap():
    result = radial_solver(...)
    k_inner = np.copy(result.k)
    return k_inner
```

Or, store the value in a stack allocated variable and return it instead (again the python interpreter will handle the
memory of the variable). This can be done by just accessing the first value of the array (if that is the only value you
need):

```python
import numpy as np

def wrap():
    result = radial_solver(...)
    k_inner = result.k[0]  # k_inner is now a double complex floating point number rather than an array.
    return k_inner
```

## Shooting Method

### Common Reasons for Integration Failure:
Below are a list of common problems that lead to integration failure. These are grouped by message codes which can be
accessed via `rs_solution.message`.

#### "Error in step size calculation:\n\tRequired step size is less than spacing between numbers."
This message indicates that the integrator could not solve the problem. 
- Possible causes:
    - `integration_rtol` or `integration_atol` is too small. Or, counterintuitively, is too large and led to compounding errors.
    - If there is a liquid layer and it is not static (via the `is_static_by_layer` variable) then the integrator likely ran into an unstable solution. This is particularly common if the planet's forcing frequency is too small (forcing periods >~ 3 days can cause this). Suggest decreasing forcing period or change the liquid layer to be static.

#### "Maximum number of steps (set by user) exceeded during integration."
This message indicates that the integrator hit the `max_num_steps` parameter during integration. You can try to increase this value but this is a good indication that the integrator is having difficulty solving the problem. Increasing the max number of steps will likely lead to a very slow integration or raise one of the other problems listed here.

#### "Maximum number of steps (set by system architecture) exceeded during integration."
This message indicates that the integrator's solution arrays exceeded the `max_ram_MB` size during integration. You can try to increase this value but this is a good indication that the integrator is having difficulty solving the problem. Increasing the max number of ram will likely lead to a very slow integration or raise one of the other problems listed here.

#### Crash
In the unlikely scenario where the integration crashes with no warnings or exceptions please try the following:
- If time permits, rerun the code to see if the crash occurs again.
    - During the rerun, keep an eye on system memory usage to ensure it is not being used at 100%.
- If crash does reoccur:
    - Ensure all arrays listed above are the same size. If any 

#### NaNs are returned
Sometimes the integration is successful but returns NaN results for all the Love Numbers. Try the following:
- Check the upper radius of each layer. NaNs may result if the specified upper radius is lower than the upper radius specified by the input radius array


## Propagation Matrix Method
An alternative to the shooting method for solving the viscoelastic-gravitational problem within planets is to use the
_Propagation Matrix_ method. This approach reduces the problem to matrix operations and does not require any direct
integration of the differential equations. The major advantage to this approach is that the number of operations is both
fixed and small, leading much faster computations. The downside is that the layer assumptions are limited. As currently
implemented in TidalPy, the propagation method does not allow for: liquid layers, compressibility, or dynamic tides. If
none of these are applicable to your problem then it is likely a much faster and stable approach. A good review of this
method is presented in [Sabadini, Vermeersen, & Cambiotti (2016)](https://www.barnesandnoble.com/w/global-dynamics-of-the-earth-roberto-sabadini/1123259823) (hereafter, SVC16).

### Solid Layer Fundamental Matrix
TidalPy implements the fundamental matrix defined by SVC16 Eq. 2.42 and its inverse Eq. 2.47 in 
`TidalPy.RadialSolver.PropMatrix.solid_matrix`. 

# References

* Numerical Shooting Method
    * [Takeuchi & Saito (1972)](https://ui.adsabs.harvard.edu/abs/1972MCPAR..11..217T/abstract)
    * [Tobie et al. (2005)](https://ui.adsabs.harvard.edu/abs/2005Icar..177..534T/abstract)
* Static / Incompressible Liquid Layers
    * [Takeuchi & Saito (1972)](https://ui.adsabs.harvard.edu/abs/1972MCPAR..11..217T/abstract)
* Dynamic / Incompressible Layers
    * [Kamata et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015JGRE..120.1528K/abstract)
    * [Beuthe (2015)](https://ui.adsabs.harvard.edu/abs/2015Icar..258..239B/abstract)
* Propagator Matrix Technique
    * [Roberts & Nimmo (2008)](https://ui.adsabs.harvard.edu/abs/2008Icar..194..675R/abstract)
    * [Henning & Hurford (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...789...30H/abstract)


## Radial Solver Solution
`TidalPy.RadialSolver.radial_solver` functions stores the viscoelastic-gravitational solution, results of solving a 
planet's equation of state, and other parameters and meta data in a cythonized python class `TidalPy.RadialSolver.rs_solution.RadialSolverSolution`. This document details how to access this API from both Python and Cython.

### Python API
The example below outlines which parameters are available via Python from the `RadialSolverSolution` class instance.

```python

rs_solution = radial_solver(...)  # See inputs from the previous section

# Check is the solution was successful.
# See the following sections for information on why a solution may not be successful.
rs_solution.success     # == (True / False)
rs_solution.error_code  # Equals 0 if no error.

# Check any integration messages (e.g., hints about why the solver failed)
rs_solution.message

# Parameters held by the equation of state solution
# See more details about these attributes in the "Equation of State Solver.md" documentation.
rs_solution.eos_error_code
rs_solution.eos_message
rs_solution.eos_success
rs_solution.eos_pressure_error
rs_solution.eos_iterations
rs_solution.eos_steps_taken

# Physical parameters held by the EOS solution
# Arrays defined radially throughout the planet
rs_solution.radius_array
rs_solution.gravity_array
rs_solution.pressure_array
rs_solution.mass_array
rs_solution.moi_array
rs_solution.density_array
rs_solution.shear_modulus_array  # Complex-valued
rs_solution.bulk_modulus_array   # Complex-valued

# Scalars (usually defined at surface or core)
rs_solution.radius
rs_solution.mass
rs_solution.moi
rs_solution.moi_factor
rs_solution.central_pressure
rs_solution.surface_pressure
rs_solution.surface_gravity

# Radial solver outputs
rs_solution.degree_l  # Which harmonic degree was used to perform calculations.

# TidalPy's radial solver uses an integration technique which solves ODEs throughout each layer of a planet. The number of integration steps is variable and determined by the local error in the calculations.
# The larger the number of steps, the lower the error but the larger the computation time. You can control what error level is acceptable by adjusting `radial_solver`s `integration_rtol`, `integration_atol`, and the `integration_method` used. Furthermore, you can force `radial_solver` to automatically fail if too many integration steps are taken by adjusting `max_num_steps` and `max_ram_MB` (useful if running MCMCs that may enter a bad parameter space leading to unstable results). 
# The solution records the number of steps taken for each solution (up to 3) in each layer. This is a np.ndarray of ints structured as [num_layers, 3]. Note that the second dimension is for the solution number. Solid layers will have 3 solutions, dynamic liquid layers will have 2, and static liquid layers will only have 1; the unused solutions will be set to 0.
# This parameter is useful to look at to both tweak tolerances to maximize performance and to check for instability (generally speaking, radial solver should not need more than ~100 steps per solution per layer) in complex scenarios this may rise to ~1000 but if you see >10_000 then it is very likely that the solution is unstable. See `rs_solution.plot_ys()` for another way to check for instability.
rs_solution.steps_taken

# All radial function results are stored in the `result` array. This is a M x 6 x N where N is the number of radial steps, M is the number of solution types (e.g., "Tidal", "Loading", etc.)
# The "6" represents the 6 different radial functions. For Solid layers these are the traditional y1, y2, y3 ... (TidalPy uses the Takeuchi & Saito (1972) format of these values).
# Liquid layers do not define y4 so it will always be nan in these layers. Static liquid layers further do not have a defined y2, y3, y5, y6. TidalPy uses the method of Saito (1974) to define a new "y7" which takes the place of "y6" in this array; all other values should be nan. 
rs_solution.result

# The radial functions can be quickly plotted as a function of radius using the following method. This is a useful debugging tool to check for instability.
# Instability generally looks like very large spikes, very constant sinusoidal results, or just results that do not smoothly change with radius.
# If multiple solution types were requested, then this will plot each on the same plots as a different color.
# This function returns a tuple containing the matplotlib (`Figure`, array of `Axis`).
rs_solution.plot_ys()

# There is a helpful way to directly get the radial solution results for a specified solution type, if multiple solutions were requested.
rs_solution['tidal']    # Tidal gravitational-viscoelastic results
rs_solution['loading']  # Loading gravitational-viscoelastic results

# Complex Love numbers can be accessed using the `love` attribute which is a np.ndarray (complex) with the shape of [num_solve_for, 3]
# The first dimension is solution type. If you are only solving for tidal Love numbers then there will only be one solution in the "0" location.
# tidal will be at 0 and load at 1). The second dimension is used to access the three Love numbers:
# - index 0: k Love number.
# - index 1: h Love number.
# - index 2: l Shida number.
# For example, if you only set `solve_for=('tidal',)` then `rs_solution.love[0, 0]` == tidal k, `rs_solution.love[0, 1]` == tidal h.
rs_solution.love

# There are also helpful shortcuts for each individual Love number.
# These are complex scalars if only one solution type is requested. Otherwise they are np.ndarrays.
rs_solution.k
rs_solution.h
rs_solution.l
```

## Helper Functions

TidalPy's radial solver function is particular on the format of its inputs. It can be easy to make a mistake which leads to solution failures, exceptions, or crashes[^1].
The helper functions in this section are designed to give users an easy interface to provide data that is then translated into the inputs required by the `radial_solver` function.

[^1]: If a crash does occur, please report it on TidalPy's GitHub issues page. Please include the exact inputs used.

### `build_planet_constant_layers`

Import with `from TidalPy.RadialSolver import build_rs_input_homogenous_layers`

Creates radial solver inputs based on user provided parameters for a planet with homogenous layers (each layer has a constant density, viscosity, shear, etc.).
Checks will be performed to ensure that the inputs are valid.

Arguments and use case:
```python
output = build_rs_input_homogenous_layers(
    planet_radius,                    # Radius of planet (float64) [m]
    forcing_frequency,                # Forcing frequency, used to solve for the complex shear / bulk modulus (float64) [rad s-1]
    density_tuple,                    # Tuple of floats for each layer's constant density. (Tuple[float64]; len = num_layers)
    static_bulk_modulus_tuple,        # Tuple of floats for each layer's constant static bulk modulus. (Tuple[float64]; len = num_layers)
    static_shear_modulus_tuple,       # Tuple of floats for each layer's constant static shear modulus. (Tuple[float64]; len = num_layers)
    bulk_viscosity_tuple,             # Tuple of floats for each layer's constant bulk viscosity. (Tuple[float64]; len = num_layers)
    shear_viscosity_tuple,            # Tuple of floats for each layer's constant shear viscosity. (Tuple[float64]; len = num_layers)
    layer_type_tuple,                 # Tuple of strings for each layer type. (Tuple[str]; len = num_layers)
    layer_is_static_tuple,            # Tuple of booleans for if each layer should use the static assumption. (Tuple[bool]; len = num_layers)
    layer_is_incompressible_tuple,    # Tuple of booleans for if each layer should use the incompressible assumption. (Tuple[bool]; len = num_layers)
    shear_rheology_model_tuple,       # Tuple of rheology instances for each layer's complex shear calculation. (Tuple[RheologyModelBase]; len = num_layers)
    bulk_rheology_model_tuple,        # Tuple of rheology instances for each layer's complex bulk calculation. (Tuple[RheologyModelBase]; len = num_layers)
    # One of the following tuples must be provided. They define the layer geometry
    radius_fraction_tuple = None,     # Tuple of floats for each layer's radius fraction (R_layer / R_Planet).  (Tuple[float64]; len = num_layers)
    thickness_fraction_tuple = None,  # Tuple of floats for each layer's thickness fraction ((R_layer - R_layer_below) / R_Planet).  (Tuple[float64]; len = num_layers)
    volume_fraction_tuple = None,     # Tuple of floats for each layer's volume fraction (V_layer / V_Planet).  (Tuple[float64]; len = num_layers)
    # Each layer is further sub-divided into slices. At least 5 slices per layer is required. 
    slices_tuple = None,              # Tuple of ints for the number of subslices each layer should be built with. (Tuple[int]; len = num_layers)
    slice_per_layer = 10,             # Number of slices to build all layers with. Used for each layer if `slices_tuple` is not provided (int)
    cpp_bool perform_checks = True    # Flag to tell function to perform additional checks on inputs. There is a small performance hit, but recommended unless you are sure your input is valid. (boolean)

# The output is a named tuple with the following attributes:
output.radius_array
output.density_array
output.complex_bulk_modulus_array
output.complex_shear_modulus_array
output.frequency
output.planet_bulk_density
output.layer_types
output.is_static_bylayer
output.is_incompressible_bylayer
output.upper_radius_bylayer_array

# These are all of the required, positional arguments of TidalPy's `radial_solver` function. 
# They can be easily passed to the solver
solution = radial_solver(
    *output
    # Any changes to keyword arguments here....
    )
```
