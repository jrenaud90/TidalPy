# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

cimport numpy as cnp
cnp.import_array()

import numpy as np

from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, make_unique
from libcpp cimport bool as cpp_bool

from CyRK cimport PreEvalFunc, ODEMethod

from TidalPy.constants cimport tidalpy_config_ptr, get_shared_config_address

# Wire up the pointer at import time
if tidalpy_config_ptr == NULL:
    tidalpy_config_ptr = get_shared_config_address()

from TidalPy.Material_x.eos.ode cimport c_EOS_ODEInput, C_EOS_Y_VALUES, C_EOS_EXTRA_VALUES, C_EOS_DY_VALUES
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.Material_x.eos.methods.interpolate cimport c_InterpolateEOSInput, c_preeval_interpolate


def solve_eos(
        tuple radius_array_bylayer,
        tuple density_array_bylayer,
        tuple complex_bulk_modulus_array_bylayer,
        tuple complex_shear_modulus_array_bylayer,
        double planet_bulk_density,
        tuple upper_radius_bylayer,
        double surface_pressure = 0.0,
        double G_to_use = -1.0,
        str integration_method = 'DOP853',
        double rtol = 1.0e-6,
        double atol = 1.0e-10,
        double pressure_tol = 1.0e-3,
        size_t max_iters = 100,
        bint verbose = True
    ):
    """
    Solve the equation of state for a layered planet.

    Parameters
    ----------
    radius_array_bylayer : tuple of ndarray[double]
        Radial positions from center to surface [m] for each layer.
    density_array_bylayer : tuple of ndarray[double]
        Density at each radius [kg m-3] for each layer.
    complex_bulk_modulus_array_bylayer : tuple of ndarray[complex128]
        Complex bulk modulus at each radius [Pa] for each layer.
    complex_shear_modulus_array_bylayer : tuple of ndarray[complex128]
        Complex shear modulus at each radius [Pa] for each layer.
    planet_bulk_density : double
        Bulk density of the planet [kg m-3].
    upper_radius_bylayer : tuple of float
        Upper radius of each macro layer [m].
    surface_pressure : double, optional
        Expected surface pressure [Pa]. Default 0.0.
    G_to_use : double, optional
        Gravitational constant [m3 kg-1 s-2]. If -1.0, uses TidalPy config value.
    integration_method : str, optional
        Integration method name. Default 'DOP853'.
    rtol : double, optional
        Relative tolerance. Default 1e-6.
    atol : double, optional
        Absolute tolerance. Default 1e-10.
    pressure_tol : double, optional
        Pressure convergence tolerance. Default 1e-3.
    max_iters : int, optional
        Maximum convergence iterations. Default 100.
    verbose : bool, optional
        Print status messages. Default True.

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'success' : bool
        - 'message' : str
        - 'iterations' : int
        - 'pressure_error' : float
        - 'gravity' : ndarray
        - 'pressure' : ndarray
        - 'mass' : ndarray
        - 'moi' : ndarray
        - 'density' : ndarray
        - 'complex_shear' : ndarray[complex128]
        - 'complex_bulk' : ndarray[complex128]
        - 'surface_gravity' : float
        - 'surface_pressure' : float
        - 'central_pressure' : float
        - 'planet_mass' : float
        - 'planet_moi' : float
    """

    # Use config G if not provided
    if G_to_use < 0.0:
        G_to_use = tidalpy_config_ptr.d_G

    # Parse integration method
    cdef ODEMethod ode_method
    cdef str method_upper = integration_method.upper()
    if method_upper == 'DOP853':
        ode_method = ODEMethod.DOP853
    elif method_upper == 'RK45':
        ode_method = ODEMethod.RK45
    elif method_upper == 'RK23':
        ode_method = ODEMethod.RK23
    else:
        raise ValueError(f"Unsupported integration method: {integration_method}")

    # Setup layer information
    cdef size_t num_layers = len(upper_radius_bylayer)
    if num_layers == 0:
        raise ValueError("Must provide at least one layer in upper_radius_bylayer.")

    # Build upper radius array
    cdef cnp.ndarray[double, ndim=1, mode='c'] upper_radius_arr = np.array(upper_radius_bylayer, dtype=np.float64)
    cdef double* upper_radius_ptr = &upper_radius_arr[0]

    # Build full radius array for the C++ EOSSolution object
    cdef cnp.ndarray[double, ndim=1, mode='c'] full_radius_array = np.concatenate(radius_array_bylayer).astype(np.float64)
    cdef size_t total_slices = full_radius_array.shape[0]
    cdef double* radius_ptr = &full_radius_array[0]

    # Create the EOS solution object
    cdef unique_ptr[c_EOSSolution] eos_solution_uptr = make_unique[c_EOSSolution](
        upper_radius_ptr,
        num_layers,
        radius_ptr,
        total_slices
    )
    cdef c_EOSSolution* eos_solution_ptr = eos_solution_uptr.get()

    # Build EOS function and input vectors for each layer (all using interpolation method)
    cdef vector[PreEvalFunc] eos_function_vec
    cdef vector[c_EOS_ODEInput] eos_input_vec
    
    # Pre-allocate a C++ vector to store the interp structs so their memory addresses stay valid safely
    cdef vector[c_InterpolateEOSInput] interp_inputs_vec
    interp_inputs_vec.resize(num_layers)

    cdef c_EOS_ODEInput eos_ode_input
    eos_ode_input.G_to_use      = G_to_use
    eos_ode_input.planet_radius = upper_radius_arr[num_layers - 1]
    eos_ode_input.final_solve   = False
    eos_ode_input.update_bulk   = False
    eos_ode_input.update_shear  = False

    # Define memoryviews to safely extract continuous pointers from the python tuples
    cdef double[::1] r_view
    cdef double[::1] d_view
    cdef double complex[::1] b_view
    cdef double complex[::1] s_view

    cdef size_t layer_i
    for layer_i in range(num_layers):
        # Extract the arrays for this specific layer
        r_view = radius_array_bylayer[layer_i]
        d_view = density_array_bylayer[layer_i]
        b_view = complex_bulk_modulus_array_bylayer[layer_i]
        s_view = complex_shear_modulus_array_bylayer[layer_i]

        interp_inputs_vec[layer_i].num_slices              = r_view.shape[0]
        interp_inputs_vec[layer_i].radius_array_ptr        = &r_view[0]
        interp_inputs_vec[layer_i].density_array_ptr       = &d_view[0]
        interp_inputs_vec[layer_i].bulk_modulus_array_ptr  = <cpp_complex[double]*>&b_view[0]
        interp_inputs_vec[layer_i].shear_modulus_array_ptr = <cpp_complex[double]*>&s_view[0]

        # Link the ODE input to this specific layer's struct memory
        eos_ode_input.eos_input_ptr = <char*>(&interp_inputs_vec[layer_i])

        eos_function_vec.push_back(c_preeval_interpolate)
        eos_input_vec.push_back(eos_ode_input)

    # Call the C++ solver
    with nogil:
        c_solve_eos(
            eos_solution_ptr,
            eos_function_vec,
            eos_input_vec,
            planet_bulk_density,
            surface_pressure,
            G_to_use,
            ode_method,
            rtol,
            atol,
            pressure_tol,
            max_iters,
            <cpp_bool>verbose
        )

    # Convert results to Python
    cdef size_t n = eos_solution_ptr.radius_array_size

    # Build output arrays
    cdef cnp.ndarray gravity_out  = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray pressure_out = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray mass_out     = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray moi_out      = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray density_out  = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray shear_out    = np.empty(n, dtype=np.complex128)
    cdef cnp.ndarray bulk_out     = np.empty(n, dtype=np.complex128)

    if eos_solution_ptr.success and eos_solution_ptr.other_vecs_set:
        # Copy data from C++ vectors to numpy arrays
        for i in range(n):
            gravity_out[i]  = eos_solution_ptr.gravity_array_vec[i]
            pressure_out[i] = eos_solution_ptr.pressure_array_vec[i]
            mass_out[i]     = eos_solution_ptr.mass_array_vec[i]
            moi_out[i]      = eos_solution_ptr.moi_array_vec[i]
            density_out[i]  = eos_solution_ptr.density_array_vec[i]
            shear_out[i]    = <double complex>eos_solution_ptr.complex_shear_array_vec[i]
            bulk_out[i]     = <double complex>eos_solution_ptr.complex_bulk_array_vec[i]

    return {
        'success'          : eos_solution_ptr.success,
        'message'          : eos_solution_ptr.message.decode('utf-8') if isinstance(eos_solution_ptr.message, bytes) else str(eos_solution_ptr.message),
        'iterations'       : eos_solution_ptr.iterations,
        'pressure_error'   : eos_solution_ptr.pressure_error,
        'gravity'          : gravity_out,
        'pressure'         : pressure_out,
        'mass'             : mass_out,
        'moi'              : moi_out,
        'density'          : density_out,
        'complex_shear'    : shear_out,
        'complex_bulk'     : bulk_out,
        'surface_gravity'  : eos_solution_ptr.surface_gravity,
        'surface_pressure' : eos_solution_ptr.surface_pressure,
        'central_pressure' : eos_solution_ptr.central_pressure,
        'planet_mass'      : eos_solution_ptr.mass,
        'planet_moi'       : eos_solution_ptr.moi,
    }
