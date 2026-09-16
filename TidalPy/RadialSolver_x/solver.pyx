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

from TidalPy.Utilities_x.logging_x.logger import log_warning
from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr, tidalpy_config_ptr, TidalPyConfig
from TidalPy.constants import ODE_METHOD_NAMES
# Make sure TidalPy Config Pointer is set.
set_tidalpy_config_ptr(get_shared_config_address())

from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver_x.rs_solution cimport RadialSolverSolution
from TidalPy.RadialSolver_x.rs_solution import check_surface_solve_conditioning
from TidalPy.Tides_x.love.love cimport c_parse_love_method_int


cdef cpp_bool _resolve_prop_matrix(str love_method) except *:
    # Map a Love-number method name (or alias) onto the two radial techniques this array API offers.
    cdef int method = c_parse_love_method_int(love_method.encode('utf-8'))
    if method == 0:
        return False
    if method == 1:
        return True
    raise ValueError(
        f"The array-based radial_solver supports love_method 'radial_solver' (aliases 'shooting', 'rs') and "
        f"'propagation_matrix' (aliases 'prop_matrix', 'pm', 'prop') only; '{love_method}' needs a built world "
        f"(LayeredWorld.solve_love_numbers) or the closed-form functions in TidalPy.Tides_x.love.")


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
        start_radius_tolerance = None,
        nondimensionalize = None,
        # Shooting method parameters
        use_kamata = None,
        integration_method = None,
        integration_rtol = None,
        integration_atol = None,
        scale_rtols_bylayer_type = None,
        max_num_steps = None,
        expected_size = None,
        max_ram_MB = None,
        double max_step = 0,
        # Love-number method
        str love_method = 'radial_solver',
        # Propagation matrix method parameters
        int core_model = 0,
        # Equation of State solver parameters
        tuple eos_method_bylayer = None,
        double surface_pressure = 0.0,
        eos_integration_method = None,
        eos_rtol = None,
        eos_atol = None,
        eos_pressure_tol = None,
        eos_max_iters = None,
        # Error and log reporting
        cpp_bool verbose = False,
        cpp_bool warnings = True,
        cpp_bool raise_on_fail = False,
        cpp_bool perform_checks = True,  # Maintained for API compatibility, but C++ handles checking unconditionally
        cpp_bool log_info = False
        ):
    """
    Solves the viscoelastic-gravitational problem for a planet comprised of solid and liquid layers.

    Every solver setting left as ``None`` takes the ``[radial_solver]`` (shooting method) or ``[eos_solver]``
    value of the TidalPy configuration (``TidalPy.config_x``), the same defaults the world-attached solves use.

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
    start_radius_tolerance : float64, optional
        Tolerance of the automatic starting radius, R * tol^(1/l). None: from the configuration.
    nondimensionalize : bool, optional
        Non-dimensionalize inputs before integration (the EOS and the shooting solve). None: from the
        configuration.
    use_kamata : bool, optional
        Use Kamata+ (2015) starting conditions. None: from the configuration.
    integration_method : str, optional
        CyRK integration method: 'RK23', 'RK45', 'DOP853', 'BDF', 'LSODA', 'Radau'. None: from the configuration.
    integration_rtol : float64, optional
        Relative integration tolerance. None: from the configuration.
    integration_atol : float64, optional
        Absolute integration tolerance. None: from the configuration.
    scale_rtols_bylayer_type : bool, optional
        Scale tolerances by layer type. None: from the configuration.
    max_num_steps : uint, optional
        Maximum integration steps. None: from the configuration.
    expected_size : uint, optional
        Expected integration steps per solution. None: from the configuration.
    max_ram_MB : uint, optional
        Maximum RAM for integrator [MB]. None: from the configuration.
    max_step : float64, default=0
        Maximum step size. 0 = auto-determine.
    love_method : str, default='radial_solver'
        Radial technique: 'radial_solver' (aliases 'shooting', 'rs') integrates the radial ODEs from the
        starting radius to the surface; 'propagation_matrix' (aliases 'prop_matrix', 'pm', 'prop') uses the
        matrix method, which is only valid for a single solid, static, incompressible layer. The analytic
        methods ('homogeneous', 'cpl', 'ctl') need a built world (LayeredWorld.solve_love_numbers) or the
        closed-form functions in TidalPy.Tides_x.love.
    core_model : int, default=0
        Core model for prop matrix method (0-4).
    eos_method_bylayer : tuple, default=None
        EOS method per layer. None = "interpolation" for all.
    surface_pressure : float64, default=0.0
        Planet surface pressure [Pa].
    eos_integration_method : str, optional
        EOS integration method: 'RK23', 'RK45', 'DOP853', 'BDF', 'LSODA', or 'Radau'. None: from the
        configuration.
    eos_rtol : float64, optional
        EOS relative tolerance. None: from the configuration.
    eos_atol : float64, optional
        EOS absolute tolerance. None: from the configuration.
    eos_pressure_tol : float64, optional
        Convergence tolerance on the surface-pressure mismatch, relative to the central-pressure scale
        (2/3) pi G rho^2 R^2 (keep it above ``eos_rtol``). None: from the configuration.
    eos_max_iters : int, optional
        Maximum central-pressure iterations. None: from the configuration.
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

    # Solver settings not given here come from the [radial_solver] and [eos_solver] sections of the TidalPy
    # configuration, held by the shared C++ config (the same values the world-attached solves start from).
    cdef TidalPyConfig* shared_config = tidalpy_config_ptr
    if start_radius_tolerance is None:
        start_radius_tolerance = shared_config.d_RADIAL_SOLVER_START_RADIUS_TOL
    if nondimensionalize is None:
        nondimensionalize = shared_config.d_RADIAL_SOLVER_NONDIMENSIONALIZE
    if use_kamata is None:
        use_kamata = shared_config.d_RADIAL_SOLVER_USE_KAMATA
    if integration_method is None:
        integration_method = ODE_METHOD_NAMES[shared_config.d_RADIAL_SOLVER_METHOD]
    if integration_rtol is None:
        integration_rtol = shared_config.d_RADIAL_SOLVER_RTOL
    if integration_atol is None:
        integration_atol = shared_config.d_RADIAL_SOLVER_ATOL
    if scale_rtols_bylayer_type is None:
        scale_rtols_bylayer_type = shared_config.d_RADIAL_SOLVER_SCALE_RTOLS
    if max_num_steps is None:
        max_num_steps = shared_config.d_RADIAL_SOLVER_MAX_NUM_STEPS
    if expected_size is None:
        expected_size = shared_config.d_RADIAL_SOLVER_EXPECTED_SIZE
    if max_ram_MB is None:
        max_ram_MB = shared_config.d_RADIAL_SOLVER_MAX_RAM_MB
    if eos_integration_method is None:
        eos_integration_method = ODE_METHOD_NAMES[shared_config.d_EOS_SOLVER_METHOD]
    if eos_rtol is None:
        eos_rtol = shared_config.d_EOS_SOLVER_RTOL
    if eos_atol is None:
        eos_atol = shared_config.d_EOS_SOLVER_ATOL
    if eos_pressure_tol is None:
        eos_pressure_tol = shared_config.d_EOS_SOLVER_PRESSURE_TOL
    if eos_max_iters is None:
        eos_max_iters = shared_config.d_EOS_SOLVER_MAX_ITERS

    cdef double   c_start_radius_tolerance = <double>start_radius_tolerance
    cdef cpp_bool c_nondimensionalize      = <cpp_bool>bool(nondimensionalize)
    cdef cpp_bool c_use_kamata             = <cpp_bool>bool(use_kamata)
    cdef str      c_integration_method     = str(integration_method)
    cdef double   c_integration_rtol       = <double>integration_rtol
    cdef double   c_integration_atol       = <double>integration_atol
    cdef cpp_bool c_scale_rtols            = <cpp_bool>bool(scale_rtols_bylayer_type)
    cdef size_t   c_max_num_steps          = <size_t>int(max_num_steps)
    cdef size_t   c_expected_size          = <size_t>int(expected_size)
    cdef size_t   c_max_ram_MB             = <size_t>int(max_ram_MB)
    cdef str      c_eos_integration_method = str(eos_integration_method)
    cdef double   c_eos_rtol               = <double>eos_rtol
    cdef double   c_eos_atol               = <double>eos_atol
    cdef double   c_eos_pressure_tol       = <double>eos_pressure_tol
    cdef int      c_eos_max_iters          = <int>int(eos_max_iters)

    cdef size_t total_slices = radius_array.shape[0]
    cdef size_t num_layers   = len(layer_types)

    # Every radial array must match the radius array's length; the C++ pipeline reads and
    # writes all of them over the radius-derived slice count, so a shorter array would be
    # accessed out of bounds.
    if total_slices == 0:
        raise ValueError('radius_array must not be empty.')
    if (<size_t>density_array.shape[0] != total_slices
            or <size_t>complex_bulk_modulus_array.shape[0] != total_slices
            or <size_t>complex_shear_modulus_array.shape[0] != total_slices):
        raise ValueError(
            'density, complex bulk modulus, and complex shear modulus arrays must all match '
            f'the radius array length ({total_slices}); got {density_array.shape[0]}, '
            f'{complex_bulk_modulus_array.shape[0]}, {complex_shear_modulus_array.shape[0]}.')
    if <size_t>upper_radius_bylayer_array.shape[0] != num_layers:
        raise ValueError(
            f'upper_radius_bylayer_array length ({upper_radius_bylayer_array.shape[0]}) must match '
            f'the number of layers ({num_layers}).')
    if len(is_static_bylayer) != num_layers or len(is_incompressible_bylayer) != num_layers:
        raise ValueError('layer_types, is_static_bylayer, and is_incompressible_bylayer must have equal lengths.')

    # Convert Python tuples into C++ std::vectors of strings
    cdef vector[cpp_string] c_layer_types
    for lt in layer_types:
        c_layer_types.push_back(lt.encode('utf-8'))

    cdef vector[cpp_string] c_solve_for
    if solve_for is not None:
        # The validator writes one boundary-condition model per entry into a fixed int[5] buffer.
        if len(solve_for) > 5:
            raise ValueError(
                f'radial_solver supports at most 5 simultaneous solve_for entries; got {len(solve_for)}.')
        for sf in solve_for:
            c_solve_for.push_back(sf.encode('utf-8'))

    cdef vector[cpp_string] c_eos_method_bylayer
    if eos_method_bylayer is not None:
        for em in eos_method_bylayer:
            c_eos_method_bylayer.push_back(em.encode('utf-8'))

    cdef cpp_bool use_prop_matrix = _resolve_prop_matrix(love_method)
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
            c_integration_method.encode('utf-8'),
            c_eos_method_bylayer,
            c_eos_integration_method.encode('utf-8'),
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
            c_use_kamata,
            starting_radius,
            c_start_radius_tolerance,
            integration_method_out,
            c_integration_rtol,
            c_integration_atol,
            c_scale_rtols,
            c_max_num_steps,
            c_expected_size,
            c_max_ram_MB,
            max_step,
            c_nondimensionalize,
            use_prop_matrix,
            eos_integration_method_int_bylayer_out.data(),
            eos_integration_method_out,
            c_eos_rtol,
            c_eos_atol,
            c_eos_pressure_tol,
            c_eos_max_iters,
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
            log_warning(
                f"Large number of steps taken found in radial solver solution "
                f"(max = {np.max(solution.steps_taken)}).")
        check_surface_solve_conditioning(solution.surface_solve_amplification, c_integration_rtol)

    return solution
