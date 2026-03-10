# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Top-level radial solver entry point.
# Orchestrates EOS → shooting/matrix → love number computation.

from libc.stdlib cimport malloc, free
from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

from CyRK cimport ODEMethod

from TidalPy.logger import get_logger
from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr
# Make sure TidalPy Config Pointer is set.
set_tidalpy_config_ptr(get_shared_config_address())

from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver_x.rs_solution cimport RadialSolverSolution

log = get_logger("TidalPy")

def radial_solver(
        double[::1] radius_array,
        double[::1] density_array,
        double complex[::1] complex_bulk_modulus_array,
        double complex[::1] complex_shear_modulus_array,
        double frequency,
        double planet_bulk_density,
        tuple layer_types,
        tuple is_static_bylayer,
        tuple is_incompressible_bylayer,
        double[::1] upper_radius_bylayer_array,
        int degree_l = 2,
        tuple solve_for = None,
        double starting_radius = 0.0,
        double start_radius_tolerance = 1.0e-5,
        cpp_bool nondimensionalize = True,
        # Shooting method parameters
        cpp_bool use_kamata = False,
        str integration_method = 'DOP853',
        double integration_rtol = 1.0e-5,
        double integration_atol = 1.0e-8,
        cpp_bool scale_rtols_bylayer_type = False,
        size_t max_num_steps = 500_000,
        size_t expected_size = 1000,
        size_t max_ram_MB = 500,
        double max_step = 0,
        # Propagation matrix method parameters
        cpp_bool use_prop_matrix = False,
        int core_model = 0,
        # Equation of State solver parameters
        tuple eos_method_bylayer = None,
        double surface_pressure = 0.0,
        str eos_integration_method = 'DOP853',
        double eos_rtol = 1.0e-3,
        double eos_atol = 1.0e-5,
        double eos_pressure_tol = 1.0e-3,
        int eos_max_iters = 40,
        # Error and log reporting
        cpp_bool verbose = False,
        cpp_bool warnings = True,
        cpp_bool raise_on_fail = False,
        cpp_bool perform_checks = True,  # Maintained for API compatibility, but C++ handles checking unconditionally
        cpp_bool log_info = False
        ):
    """
    Solves the viscoelastic-gravitational problem for a planet comprised of solid and liquid layers.

    Parameters
    ----------
    radius_array : np.ndarray[dtype=np.float64]
        Radius values defined at slices throughout the planet [m].
    density_array : np.ndarray[dtype=np.float64]
        Density at each radius [kg m-3].
    complex_bulk_modulus_array : np.ndarray[dtype=np.complex128]
        Bulk modulus at each radius [Pa].
    complex_shear_modulus_array : np.ndarray[dtype=np.complex128]
        Complex shear modulus at each radius [Pa].
    frequency : float64
        Forcing frequency [rad s-1]
    planet_bulk_density : float64
        Bulk density of the planet [kg m-3].
    layer_types : tuple[string, ...]
        Indicator of layer type: "solid" or "liquid".
    is_static_bylayer : tuple[bool, ...]
        Flag declaring if each layer uses the static (True) or dynamic (False) assumption.
    is_incompressible_bylayer : tuple[bool, ...]
        Flag declaring if each layer is incompressible (True) or compressible (False).
    upper_radius_bylayer_array : np.ndarray[dtype=np.float64]
        Upper radius of each layer.
    degree_l : int, default=2
        Harmonic degree.
    solve_for : tuple[str, ...], default=None
        Tuple of requested solutions ("tidal", "loading", "free"). None defaults to ("tidal",).
    starting_radius : float64, default=0.0
        Starting radius [m]. 0.0 = auto-determine.
    start_radius_tolerance : float64, default=1.0e-5
        Tolerance for starting radius formula.
    nondimensionalize : bool, default=True
        Non-dimensionalize inputs before integration.
    use_kamata : bool, default=False
        Use Kamata+ (2015) starting conditions.
    integration_method : str, default='DOP853'
        CyRK integration method ('RK23', 'RK45', 'DOP853').
    integration_rtol : float64, default=1.0e-5
        Relative integration tolerance.
    integration_atol : float64, default=1.0e-8
        Absolute integration tolerance.
    scale_rtols_bylayer_type : bool, default=False
        Scale tolerances by layer type.
    max_num_steps : uint, default=500000
        Maximum integration steps.
    expected_size : uint, default=1000
        Expected integration steps per solution.
    max_ram_MB : uint, default=500
        Maximum RAM for integrator [MB].
    max_step : float64, default=0
        Maximum step size. 0 = auto-determine.
    use_prop_matrix : bool, default=False
        Use propagation matrix method.
    core_model : int, default=0
        Core model for prop matrix method (0-4).
    eos_method_bylayer : tuple, default=None
        EOS method per layer. None = "interpolation" for all.
    surface_pressure : float64, default=0.0
        Planet surface pressure [Pa].
    eos_integration_method : str, default='DOP853'
        EOS integration method.
    eos_rtol : float64, default=1.0e-3
        EOS relative tolerance.
    eos_atol : float64, default=1.0e-5
        EOS absolute tolerance.
    eos_pressure_tol : float64, default=1.0e-3
        Pressure convergence tolerance.
    eos_max_iters : int, default=40
        Maximum EOS convergence iterations.
    verbose : bool, default=False
        Print status messages.
    warnings : bool, default=True
        Print warnings.
    raise_on_fail : bool, default=False
        Raise exception on failure.
    perform_checks : bool, default=True
        Perform input sanity checks.
    log_info : bool, default=False
        Log diagnostic information.

    Returns
    -------
    solution : RadialSolverSolution
    """

    cdef size_t total_slices = radius_array.shape[0]
    cdef size_t num_layers   = len(layer_types)

    # Convert Python tuples into C++ std::vectors of strings
    cdef vector[cpp_string] c_layer_types
    for lt in layer_types:
        c_layer_types.push_back(lt.encode('utf-8'))

    cdef vector[cpp_string] c_solve_for
    if solve_for is not None:
        for sf in solve_for:
            c_solve_for.push_back(sf.encode('utf-8'))

    cdef vector[cpp_string] c_eos_method_bylayer
    if eos_method_bylayer is not None:
        for em in eos_method_bylayer:
            c_eos_method_bylayer.push_back(em.encode('utf-8'))

    cdef vector[int] layer_types_out = vector[int](num_layers)
    # Use standard malloc for boolean arrays since std::vector<bool> behaves like a bitfield in C++
    cdef cpp_bool* c_is_static = <cpp_bool*>malloc(num_layers * sizeof(cpp_bool))
    cdef cpp_bool* c_is_incomp = <cpp_bool*>malloc(num_layers * sizeof(cpp_bool))
    
    if not c_is_static or not c_is_incomp:
        raise MemoryError("Failed to allocate memory for boolean assumption arrays.")

    for i in range(num_layers):
        if is_static_bylayer[i]:
            c_is_static[i] = True
        else:
            c_is_static[i] = False
        if is_incompressible_bylayer[i]:
            c_is_incomp[i] = True
        else:
            c_is_incomp[i] = False

    # Prepare structures for C++ validator outputs
    cdef int[5] bc_models_out
    cdef size_t num_bc_models_out = 0
    cdef ODEMethod integration_method_out     = ODEMethod.NO_METHOD_SET
    cdef ODEMethod eos_integration_method_out = ODEMethod.NO_METHOD_SET
    cdef vector[int] eos_integration_method_int_bylayer_out
    
    cdef int rs_error_code = 0
    cdef RadialSolverSolution solution
    
    try:
        # Call C++ helper to validate data and prep variables (throws ValueError on checks failure)
        c_validate_and_prep_radial_inputs(
            total_slices,
            &radius_array[0],
            &density_array[0],
            frequency,
            num_layers,
            c_layer_types,
            c_is_static,
            c_is_incomp,
            &upper_radius_bylayer_array[0],
            use_prop_matrix,
            starting_radius,
            c_solve_for,
            integration_method.encode('utf-8'),
            c_eos_method_bylayer,
            eos_integration_method.encode('utf-8'),
            warnings,
            layer_types_out.data(),
            &bc_models_out[0],
            num_bc_models_out,
            integration_method_out,
            eos_integration_method_int_bylayer_out,
            eos_integration_method_out
        )

        # Build solution storage
        solution = RadialSolverSolution(
            num_bc_models_out,
            upper_radius_bylayer_array,
            radius_array,
            degree_l
        )

        solution.set_model_names(&bc_models_out[0])

        # Run C++ radial solver 
        rs_error_code = c_radial_solver(
            solution.solution_storage_uptr.get(),
            total_slices,
            &radius_array[0],
            &density_array[0],
            <cpp_complex[double]*>&complex_bulk_modulus_array[0],
            <cpp_complex[double]*>&complex_shear_modulus_array[0],
            frequency,
            planet_bulk_density,
            num_layers,
            layer_types_out.data(),
            c_is_static,
            c_is_incomp,
            surface_pressure,
            degree_l,
            num_bc_models_out,
            &bc_models_out[0],
            core_model,
            use_kamata,
            starting_radius,
            start_radius_tolerance,
            integration_method_out,
            integration_rtol,
            integration_atol,
            scale_rtols_bylayer_type,
            max_num_steps,
            expected_size,
            max_ram_MB,
            max_step,
            nondimensionalize,
            use_prop_matrix,
            eos_integration_method_int_bylayer_out.data(),
            eos_integration_method_out,
            eos_rtol,
            eos_atol,
            eos_pressure_tol,
            eos_max_iters,
            verbose,
            warnings
        )
        
    finally:
        # Guarantee cleanup of manually allocated arrays
        free(c_is_static)
        free(c_is_incomp)

    # Finalize
    solution.finalize_python_storage()

    if log_info:
        solution.print_diagnostics(print_diagnostics=False, log_diagnostics=True)

    if ((not solution.success) or (rs_error_code < 0)) and raise_on_fail:
        if "not implemented" in solution.message.lower():
            raise NotImplementedError(solution.message)
        else:
            raise SolutionFailedError(solution.message)

    if warnings:
        if np.any(solution.steps_taken > 7_000):
            log.warning(f"Large number of steps taken found in radial solver solution (max = {np.max(solution.steps_taken)}).")

    return solution
