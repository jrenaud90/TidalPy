# Calculating Love Numbers using TidalPy's RadialSolver module and functions
_Instructions and examples on how to use TidalPy's `RadialSolver` module._

TidalPy's `RadialSolver` package allows a user to estiamte a planet's global, viscoelastic
[Love numbers](https://en.wikipedia.org/wiki/Love_number). These numbers can then be used to determine the magnitude of tidal
dissipation, speed of rotational/orbital changes, and provide predictions for gravity field measurements. 

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
radial_solver_solution = radial_solver(
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
accessed via `radial_solver_solution.message`.

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