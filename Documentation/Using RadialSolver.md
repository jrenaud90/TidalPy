# Radial Solver Documents
_Instructions and examples on how to use TidalPy's `RadialSolver` module._

## `radial_solver` Function
The `radial_solver` function, contained in the `TidalPy.RadialSolver` module can be used with the following arguments.

```python
from TidalPy.RadialSolver import radial_solver
radial_solver_solution = radial_solver(
    # Following arrays should define the respective physical properties throughout the planet. 
    # It is recommended to use at least 20 slices per layer (so the arrays should be 20+ x num_layers in size).
    # These arrays can be numpy np.ndarrays.
    radius_array,  # [m] (type: double array)
    density_array,  # [kg m-3] (type: double array)
    gravity_array,  # [m s-2] (type: double array)
    bulk_modulus_array,  # [Pa] (type: double array)
    # Note that `complex_shear_modulus_array` is the only complex-valued array. It is the result of applying a 
    #  rheological function to the planet's shear modulus and viscosity.
    complex_shear_modulus_array,  # [Pa] (double complex array)
    frequency,  # Scalar forcing frequency [rad s-1] (type: double)
    planet_bulk_density,  # Scalar bulk density of planet [kg m-3] (type: double)
    # The following tuples define the assumptions used for the major layers within a planet.
    is_solid_by_layer,  # Is each layer solid (type: tuple of bools)
    is_static_by_layer,  # Is each layer using the static tidal assumption (type: tuple of bools)
    is_incompressible_by_layer,  # Is each layer using the incompressible tidal assumption (type: tuple of bools)
    upper_radius_by_layer,  # Upper radius of each layer [m] (type: tuple of doubles)
    degree_l = 2,  # Harmonic degree used
    tuple solve_for = None,
    bool_cpp_t use_kamata = False,
    int integration_method = 1,
    double integration_rtol = 1.0e-4,
    double integration_atol = 1.0e-12,
    bool_cpp_t scale_rtols_by_layer_type = True,
    size_t max_num_steps = 500_000,
    size_t expected_size = 500,
    size_t max_ram_MB = 500,
    double max_step = 0,
    bool_cpp_t limit_solution_to_radius = True,
    bool_cpp_t nondimensionalize = True,
    bool_cpp_t verbose = False,
    bool_cpp_t raise_on_fail = False
    )

```

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
## `RadialSolverSolution` Class

### Love Numbers
In addition to the full `y` radial solution results, you can also quickly access the Love numbers using:

```python
solution.love  
# Returns all Love numbers for each solution type specified.
```

These are returned in a 2D array structured like [num_solution_types : love_number] where the first column is the
solution type. If you are only solving for tidal Love numbers then there will only be one solution in the "0" location.
If you are solving for multiple solutions then they will be stored in order (e.g., if you solve for tidal and load then
tidal will be at 0 and load at 1). The second column is used to access the three Love numbers:
- index 0: k Love number.
- index 1: h Love number.
- index 2: l Shida number.

Example:

Say you want to access the h Love number for solution 0:

```python
k = solution.love[0, 1]
```

There are shortcuts to make accessing specific Love numbers easier. The following attributes will return arrays for the
respective Love number for each solution type.

```python
solution.k = # [Solution 0's k, Solution 1's k, ...]
solution.h = # [Solution 0's h, Solution 1's h, ...]
solution.l = # [Solution 0's l, Solution 1's l, ...]
```

**Performance Note**

For all of these attributes (`.love; .k; .h; .l`) the arrays are being created each time the attribute is accessed.
If you have a code that accesses these numbers often (more than once) it is better to store them in a local variable.

```python
k = solution.k  # Creates an array for all k Love numbers.

# ... Do a bunch of stuff with the new "k" variable.
```
