
# Radial Solver Solution: `RadialSolverSolution` Class
`TidalPy.RadialSolver.radial_solver` functions stores the viscoelastic-gravitational solution, results of solving a 
planet's equation of state, and other parameters and meta data in a cythonized python class `TidalPy.RadialSolver.rs_solution.RadialSolverSolution`. This document details how to access this API from both Python and Cython.

## Python API
The example below outlines which parameters are available via Python from the `RadialSolverSolution` class instance.

```python

rs_solution = radial_solver(...)
# `radial_solver` and its arguments are discussed in "1 - Calculating Love Numbers with TidalPy.md"

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

# Layer specific data
rs_solution.layer_upper_radius_array

# Scalars (usually defined at surface or core)
rs_solution.radius
rs_solution.volume
rs_solution.mass
rs_solution.moi
rs_solution.moi_factor
rs_solution.density_bulk
rs_solution.central_pressure
rs_solution.surface_pressure
rs_solution.surface_gravity

# Radial solver outputs
rs_solution.num_ytypes  # Number of y-solutions solved simultaneously (e.g., 2 for ('tidal', 'loading'))
rs_solution.num_layers  # Number of layers
rs_solution.degree_l    # Which harmonic degree was used to perform calculations.

# All radial function results are stored in the `result` array. This is a M x 6 x N where N is the number of radial steps,
# M is the number of solution types (e.g., "Tidal", "Loading", etc.)
# The "6" represents the 6 different radial functions. For Solid layers these are the traditional y1, y2, y3 ...
# (TidalPy uses the Takeuchi & Saito (1972) format of these values).
# Liquid layers do not define y4 so it will always be nan in these layers. Static liquid layers further do not
# have a defined y2, y3, y5, y6. TidalPy uses the method of Saito (1974) to define a new "y7" which takes the place of
# "y6" in this array; all other values should be nan. 
rs_solution.result

# There is a helpful way to directly get the radial solution results for a specified solution type, if multiple solutions were requested.
rs_solution['tidal']    # Tidal gravitational-viscoelastic results
rs_solution['loading']  # Loading gravitational-viscoelastic results

# Complex Love numbers can be accessed using the `love` attribute which is a np.ndarray (complex)
# with the shape of [num_solve_for, 3]
# The first dimension is solution type. If you are only solving for tidal Love numbers then there will only be one
# solution in the "0" location. Tidal will be at 0 and Loading at 1). The second dimension is used to access the three Love numbers:
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

# TidalPy's radial solver uses an integration technique which solves ODEs throughout each layer of a planet.
# The number of integration steps is variable and determined by the local error in the calculations.
# The larger the number of steps, the lower the error but the larger the computation time.
# You can control what error level is acceptable by adjusting `radial_solver`s `integration_rtol`, `integration_atol`,
# and the `integration_method` used. Furthermore, you can force `radial_solver` to automatically fail if too many integration
# steps are taken by adjusting `max_num_steps` and `max_ram_MB` (useful if running MCMCs that may enter a bad parameter
# space leading to unstable results).
# The solution records the number of steps taken for each solution (up to 3) in each layer.
# This is a np.ndarray of ints structured as [num_layers, 3]. Note that the second dimension is for the solution number.
# Solid layers will have 3 solutions, dynamic liquid layers will have 2, and static liquid layers will only have 1;
# the unused solutions will be set to 0.
# This parameter is useful to look at to both tweak tolerances to maximize performance and to check for instability
# (generally speaking, radial solver should not need more than ~100 steps per solution per layer) in complex scenarios
# this may rise to ~1000 but if you see >10_000 then it is very likely that the solution is unstable.
# See `rs_solution.plot_ys()` for another way to check for instability.
rs_solution.steps_taken

# The radial functions can be quickly plotted as a function of radius using the following method. This is a useful debugging tool to check for instability.
# Instability generally looks like very large spikes, very constant sinusoidal results, or just results that do not smoothly change with radius.
# If multiple solution types were requested, then this will plot each on the same plots as a different color.
# This function returns a tuple containing the matplotlib (`Figure`, array of `Axis`).
rs_solution.plot_ys()  # Note () - this is a function

# There are helper functions to quickly plot the interior structure of the planet found by the EOS solver
rs_solution.plot_interior()

# The solution has a built in tool to create a string of information and data.
# This will print out top line results and help diagnose problems.
#   If `print_diagnostics` is True, then the string will be printed out to the python consol.
#   If `log_diagnostics` is True, then the string will be sent to TidalPy's log. Depending on the log settings this can be printed to the consol or saved to a file.
#      You can turn on this save to file feature by first calling `import TidalPy; TidalPy.log_to_file()`
rs_solution.print_diagnostics(print_diagnostics = True, log_diagnostics = False)

# The following method is provided so that the user can utilize the equation of state result that the EOS solver found even after the radial solver is finished.
eos_result_array = rs_solution.eos_call(radius=1.5e6)
# In the example above we want the equation of state results at the radius value of 1.5e6 [m]
# The `eos_result_array` is a np.ndarray with 9 components that correspond to the following evaluated at this radius.
# eos_result_array == [gravity, pressure, mass, moi component, density,
#                      shear modulus (real), shear modulus (imag),
#                      bulk modulus (real), bulk modulus (imag)]        
```

**Performance Note**

There is a minor overhead when accessing any of the Solver's attributes (e.g., `.result; .love; .k; .h; .l`).
If you have a code that accesses these numbers often (more than once) it is better to store them in a local variable.

```python
k_local = solution.k  # Performs a background lookup and np.ndarray operation to produce an array for all k Love numbers.

# ... Do a bunch of stuff with the new "k_local" variable.
```

## Cython API

The radial solver class can have much of its data accessed or modified via Cython. Below is a list of available attributes and methods.

```cython
rs_solution = radial_solver(...)

size_t radius_array_size
size_t num_ytypes
size_t num_layers
cpp_bool ytype_names_set
char* ytypes[5]

# Main storage container
shared_ptr[RadialSolutionStorageCC] solution_storage_sptr
RadialSolutionStorageCC* solution_storage_ptr

# Result pointers and data
cnp.ndarray full_solution_arr

# Love number information
cnp.ndarray complex_love_arr

# EOS solution arrays
cnp.ndarray radius_array_cnp
cnp.ndarray gravity_array_cnp
cnp.ndarray pressure_array_cnp
cnp.ndarray mass_array_cnp
cnp.ndarray moi_array_cnp
cnp.ndarray density_array_cnp
cnp.ndarray shear_modulus_array_cnp
cnp.ndarray bulk_modulus_array_cnp

# Shooting method diagnostics
cnp.ndarray shooting_method_steps_taken_array
cnp.ndarray eos_steps_taken_array

void finalize_python_storage()
void set_model_names(int* bc_models_ptr)
void change_radius_array(double* new_radius_array_ptr, size_t new_size_radius_array, cpp_bool array_changed)
```

The majority of the data is stored in the C++ class `RadialSolutionStorageCC` which has the following attributes available via Cython.

```cython
RadialSolutionStorageCC* rs_solution_cc = rs_solution.solution_storage_ptr

cpp_bool success
int error_code
int degree_l
char* message_ptr
size_t num_ytypes
size_t num_slices
size_t num_layers
size_t total_size
unique_ptr[EOSSolutionCC] eos_solution_uptr
vector[double] full_solution_vec
vector[double] complex_love_vec
vector[size_t] shooting_method_steps_taken_vec

EOSSolutionCC* get_eos_solution_ptr()
void change_radius_array(double* new_radius_array_ptr, size_t new_size_radius_array, cpp_bool array_changed)
void set_message(const char* new_message)
void find_love()
void dimensionalize_data(NonDimensionalScalesCC* nondim_scales, cpp_bool redimensionalize)
```
