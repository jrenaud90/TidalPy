# Calculating Love Numbers using TidalPy's RadialSolver module and functions
_Instructions and examples on how to use TidalPy's `RadialSolver` module._

TidalPy's `RadialSolver` package allows a user to estimate a planet's global, viscoelastic
[Love numbers](https://en.wikipedia.org/wiki/Love_number). These numbers can then be used to determine the magnitude of tidal
dissipation, speed of rotational/orbital changes, and provide predictions for gravity or distance measurements. 

TidalPy uses a numerical shooting methods which integrates a set of 3x6 viscoelastic-gravitational ordinary differential equations
from the core of the planet to its surface. The final result is a super position of the three different solutions where each
solutions' coefficients are determined by boundary conditions at the surface of the planet. 

To learn more about this method please review the literature cited in the references section.

## Radial Solver Function `TidalPy.RadialSolver.radial_solver`
The `radial_solver` function, contained in the `TidalPy.RadialSolver` module can be used with the following arguments.

Note that all arrays must be [C-contiguous](https://stackoverflow.com/questions/26998223/what-is-the-difference-between-contiguous-and-non-contiguous-arrays).
If you suspect that an array may not be C-contiguous you can use the numpy function `arr = np.ascontiguousarray(arr)` to ensure that they are before being passed to the rheology methods.

```python
from TidalPy.RadialSolver import radial_solver
rs_solution = radial_solver(
    # Following arrays should define the respective physical properties throughout the planet. 
    # It is recommended to use at least 20 slices per layer (so the arrays should be 20+ x num_layers in size).
    # These arrays can be numpy np.ndarrays.
    
    radius_array,
    # [m] (type: double array)
    
    density_array,
    # [kg m-3] (type: double array)
   
    gravity_array, 
    # [m s-2] (type: double array)
    
    bulk_modulus_array,
    # [Pa] (type: double array)
   
    complex_shear_modulus_array, 
    # [Pa] (double complex array)
    # Note that `complex_shear_modulus_array` is the only complex-valued array. It is the result of applying a 
    #  rheological function to the planet's shear modulus and viscosity.
    
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
    ## Must be provided but ignored when `use_propagation_matrix = True`
    
    is_static_by_layer,
    # Is each layer using the static tidal assumption (type: tuple of bools)
    ## Must be provided but ignored when `use_propagation_matrix = True`
    
    is_incompressible_by_layer,
    # Is each layer using the incompressible tidal assumption (type: tuple of bools)
    ## Must be provided but ignored when `use_propagation_matrix = True`
    
    upper_radius_by_layer,
    # Upper radius of each layer [m] (type: tuple of doubles)
    
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
    # You must provide these as a tuple of strings, even if you are only solving for one thing at a time. 
    # For example: solve_for = ('tidal',)  <-- Note the required comma.
    # It is more computationally efficient to solve for multiple Love number types at same time if you need more.
    # For example: solve_for = ('tidal', 'loading')
    # The number of items passed to this variable sets the size of `num_solve_for` discussed later in this document.

    use_propagation_matrix = False,
    # Setting this to True changes the underlying mechanism that TidalPy will use to solve for the radial solutions.
    # See `Propagation Matrix Method` (when this is set to True) and `Shooting Method` (when this is set to False; the default)
    # sections in this document for more details.
    
    # The following indented options are only used when `use_propagation_matrix=False` (Shooting Method)
        use_kamata = False,
        # Use starting conditions as defined in Kamata et al. (2015; JGR:P)
        # These are particularly stable for incompressible layers.
        # If set to False then the starting conditions from Takeuchi & Saito (1972)
        
        integration_method = 'RK45',
        # Integration method to use.
        # Options:
        #  - 'RK23'    Explicit Runge-Kutta method of order 3(2)
        #  - 'RK45'    Explicit Runge-Kutta method of order 5(4)
        #  - 'DOP853'  Explicit Runge-Kutta method of order 8
        
        integration_rtol = 1.0e-6,
        # Integration relative tolerance
        
        integration_atol = 1.0e-12,
        # Integration absolute tolerance (tolerance when y is near 0)
        
        scale_rtols_by_layer_type = False,
        # If True, then relative tolerance will be modified by the layer type (generally liquid layers require smaller rtol)
        
        max_num_steps = 500_000,
        # Maximum number of steps allowed for each integration (note that multiple integrations can occur depending on the
        # number and type of layers).
        
        expected_size = 500,
        # Expected number of steps needed to complete adaptive radial integration. It is better to overshoot this.
        
        max_ram_MB = 500,
        # Maximum amount of ram the integrator is allowed. Note that real ram usage will be larger than this.
        
        max_step = 0,
        # Maximum step size allowed. If 0 then the integrator will attempt to find a reasonable maximum step size.
        
        limit_solution_to_radius = True,
        # If True, then radial solutions are limited to the same radial steps defined in the user-provided `radius_array`.
    
    nondimensionalize = True,
    # If True, then the inputs will be non-dimensionalized before integration. Radial solutions will be
    # re-dimensionalized before solution is given back to the user.
    # Generally setting this to True will provide more solution stability.
    
    verbose = False,
    # If True, then additional information will be printed to the console at a small performance hit.
    
    raise_on_fail = False
    # If True, then any exceptions raised during integration will stop integration. Otherwise the integration will
    # complete without error but the solution will not be successful (you can check by using `result.success`).
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
