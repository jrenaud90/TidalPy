# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Standalone radial solver entry point: EOS, shooting or matrix solve, Love numbers.

from libc.stdlib cimport malloc, free
from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.string cimport string as cpp_string
from libcpp.utility cimport move
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

from CyRK cimport ODEMethod

from TidalPy.Utilities_x.logging_x.logger import log_warning
from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr, tidalpy_config_ptr, TidalPyConfig
from TidalPy.constants import ODE_METHOD_NAMES
set_tidalpy_config_ptr(get_shared_config_address())

from TidalPy.exceptions import SolutionFailedError
from TidalPy.RadialSolver_x.rs_solution cimport RadialSolverSolution, c_RadialSolutionStorage
from TidalPy.RadialSolver_x.rs_solution cimport cy_check_surface_solve_conditioning
from TidalPy.Tides_x.love.love cimport c_parse_love_method_int
# The world types and the C++ profile builder come from this module's own .pxd, which redeclares them rather
# than cimporting structures_x.worlds.layered; see the note there for why that import cannot be used.


# Nothing a caller can see reads this; it is here so a C++ message about that world says where it came from.
DEF PROFILE_WORLD_NAME = b"radial_solver_profile"


cdef class _ProfileWorldAnchor:
    """Owns the C++ world a supplied profile was built into, for as long as its solution is alive.

    The standalone API hands `radial_solver` arrays rather than a world, so the world it builds never reaches
    Python. It must outlive the solution, though, because the material provider the Love solve installed on
    the solution reads this world's solved EOS; `RadialSolverSolution._adopt` keeps this object alive for it.
    """
    cdef shared_ptr[c_LayeredWorld] world_sptr


cdef cpp_bool cy_resolve_prop_matrix(str love_method) except *:
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
        cpp_bool perform_checks = True,  # kept for API compatibility; C++ always validates
        cpp_bool log_info = False
        ):
    """
    Solve the viscoelastic-gravitational problem for a planet of solid and liquid layers.

    Every solver setting left as ``None`` takes the ``[radial_solver]`` or ``[eos_solver]`` value of
    ``TidalPy.config_x``, the same defaults the world-attached solves use.

    Parameters
    ----------
    radius_array : np.ndarray[dtype=np.float64]
        Radius at each slice [m]; interface radii appear twice.
    density_array : np.ndarray[dtype=np.float64]
        Density at each radius [kg m-3].
    complex_bulk_modulus_array, complex_shear_modulus_array : np.ndarray[dtype=np.complex128]
        Complex moduli at each radius [Pa].
    frequency : float64
        Forcing frequency [rad s-1].
    planet_bulk_density : float64
        [kg m-3].
    layer_types : tuple[str, ...]
        "solid" or "liquid" per layer.
    is_static_bylayer, is_incompressible_bylayer : tuple[bool, ...]
        Static (True) or dynamic, incompressible (True) or compressible, per layer.
    upper_radius_bylayer_array : np.ndarray[dtype=np.float64]
        Upper radius of each layer [m].
    degree_l : int, default=2
        Harmonic degree.
    solve_for : tuple[str, ...], optional
        Up to 5 of "tidal", "loading", "free"; None means ("tidal",).
    starting_radius : float64, default=0.0
        Starting radius [m]; 0 selects R * tol^(1/l) with ``start_radius_tolerance``.
    nondimensionalize : bool, optional
        Non-dimensionalize the EOS and shooting solves.
    use_kamata : bool, optional
        Use the Kamata et al. (2015) starting conditions.
    integration_method, eos_integration_method : str, optional
        CyRK method: 'RK23', 'RK45', 'DOP853', 'BDF', 'LSODA', or 'Radau'.
    integration_rtol, integration_atol, eos_rtol, eos_atol : float64, optional
        Integration tolerances.
    scale_rtols_bylayer_type : bool, optional
        Tighten the rtols of the unstable ys per layer type.
    max_num_steps, expected_size, max_ram_MB : uint, optional
        Integrator limits (RAM in MB).
    max_step : float64, default=0
        Maximum step size [m]; 0 selects one third of each layer thickness.
    love_method : str, default='radial_solver'
        'radial_solver' (aliases 'shooting', 'rs') or 'propagation_matrix' (aliases 'prop_matrix', 'pm',
        'prop'; a single solid, static, incompressible layer only). The analytic methods need a built world
        or the closed-form functions in TidalPy.Tides_x.love.
    core_model : int, default=0
        Propagation matrix core starting condition (0 to 4).
    eos_method_bylayer : tuple, optional
        Only "interpolate" is supported; None applies it to every layer.
    surface_pressure : float64, default=0.0
        [Pa].
    eos_pressure_tol : float64, optional
        Surface-pressure convergence tolerance relative to (2/3) pi G rho^2 R^2; keep it above ``eos_rtol``.
    eos_max_iters : int, optional
        Maximum central-pressure iterations.
    verbose, warnings, raise_on_fail, perform_checks, log_info : bool
        Reporting switches; ``perform_checks`` is accepted for compatibility and inputs are always validated.

    Returns
    -------
    solution : RadialSolverSolution
    """

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

    # The C++ pipeline indexes every array over the radius-derived slice count.
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

    cdef cpp_bool use_prop_matrix = cy_resolve_prop_matrix(love_method)
    cdef vector[int] layer_types_out = vector[int](num_layers)
    # malloc because std::vector<bool> is bit-packed.
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

    cdef int[5] bc_models_out
    cdef size_t num_bc_models_out = 0
    cdef ODEMethod integration_method_out     = ODEMethod.NO_METHOD_SET
    cdef ODEMethod eos_integration_method_out = ODEMethod.NO_METHOD_SET
    cdef vector[int] eos_integration_method_int_bylayer_out
    
    cdef int rs_error_code = 0
    cdef RadialSolverSolution solution

    # The world the supplied profile is built into, and the two solve configs. All C++: the profile never
    # becomes a Python object on its way to the solver.
    cdef vector[double] shear_static
    cdef vector[double] bulk_static
    cdef shared_ptr[c_LayeredWorld] world_sptr
    cdef c_LayeredWorld* world_ptr = NULL
    cdef c_WorldEOSSolveConfig eos_cfg
    cdef c_LoveSolveConfig love_cfg
    cdef unique_ptr[c_RadialSolutionStorage] storage_uptr
    cdef _ProfileWorldAnchor world_anchor
    cdef cpp_complex[double]* shear_ptr = NULL
    cdef cpp_complex[double]* bulk_ptr  = NULL
    cdef size_t slice_i
    
    try:
        # Raises ValueError on a failed check.
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

        # The supplied arrays describe a planet, so build that planet and solve it the way a built world is
        # solved: one code path for both APIs, with an interpolated material per layer and the complex moduli
        # handed in rather than derived from a rheology. The EOS interpolates the static moduli, the real
        # parts of the supplied complex ones. Copied into C++ vectors rather than NumPy views so no part of
        # the profile becomes a Python object on its way to the solver.
        shear_static.resize(total_slices)
        bulk_static.resize(total_slices)
        for slice_i in range(total_slices):
            shear_static[slice_i] = complex_shear_modulus_array[slice_i].real
            bulk_static[slice_i]  = complex_bulk_modulus_array[slice_i].real
        world_sptr = c_build_world_from_layered_profile(
            &radius_array[0],
            &density_array[0],
            shear_static.data(),
            bulk_static.data(),
            total_slices,
            &upper_radius_bylayer_array[0],
            layer_types_out.data(),
            c_is_static,
            c_is_incomp,
            num_layers,
            planet_bulk_density,
            PROFILE_WORLD_NAME)
    finally:
        free(c_is_static)
        free(c_is_incomp)

    world_ptr = world_sptr.get()

    # The config starts from the [eos_solver] section, exactly as the world-attached path does, with this
    # call's already-resolved settings written over it.
    eos_cfg = world_ptr.make_eos_solve_config()
    eos_cfg.surface_pressure   = surface_pressure
    eos_cfg.slices_per_layer   = <size_t>max(<int>(total_slices // num_layers), 5)
    eos_cfg.integration_method = eos_integration_method_out
    eos_cfg.rtol               = c_eos_rtol
    eos_cfg.atol               = c_eos_atol
    eos_cfg.pressure_tol       = c_eos_pressure_tol
    eos_cfg.max_iters          = <size_t>c_eos_max_iters
    eos_cfg.nondimensionalize  = c_nondimensionalize
    with nogil:
        world_ptr.solve_eos(eos_cfg)
    if world_ptr.get_eos_max_iters_hit():
        log_warning(
            f"The supplied profile's EOS solve stopped at max_iters = {c_eos_max_iters} with a surface-pressure "
            f"mismatch above eos_pressure_tol = {c_eos_pressure_tol:0.1e}; the profile is from the last "
            f"iteration. Raise eos_pressure_tol above the integration eos_rtol ({c_eos_rtol:0.1e}) or tighten "
            f"eos_rtol.")

    # From the supplied complex moduli rather than a layer rheology. The boundary-condition models and both
    # integration methods were resolved by the input check above, so nothing is re-parsed here.
    love_cfg.frequency          = frequency
    love_cfg.degree_l           = degree_l
    love_cfg.set_bc_models(&bc_models_out[0], num_bc_models_out)
    love_cfg.love_method        = <int>use_prop_matrix
    love_cfg.core_model         = core_model
    love_cfg.use_kamata         = c_use_kamata
    love_cfg.nondimensionalize  = c_nondimensionalize
    love_cfg.starting_radius    = starting_radius
    love_cfg.start_radius_tol   = c_start_radius_tolerance
    love_cfg.integration_method = integration_method_out
    love_cfg.rtol               = c_integration_rtol
    love_cfg.atol               = c_integration_atol
    love_cfg.scale_rtols        = c_scale_rtols
    love_cfg.max_num_steps      = c_max_num_steps
    love_cfg.expected_size      = c_expected_size
    love_cfg.max_ram_MB         = c_max_ram_MB
    love_cfg.max_step           = max_step
    love_cfg.verbose            = verbose
    # This function runs its own conditioning check on the finished solution below.
    love_cfg.warnings           = False

    shear_ptr = <cpp_complex[double]*><void*>&complex_shear_modulus_array[0]
    bulk_ptr  = <cpp_complex[double]*><void*>&complex_bulk_modulus_array[0]
    with nogil:
        world_ptr.solve_love_numbers_supplied(
            love_cfg,
            shear_ptr,
            bulk_ptr,
            &radius_array[0],
            total_slices)

    storage_uptr = world_ptr.release_radial_storage()
    if not storage_uptr:
        raise SolutionFailedError(
            "The radial solve produced no solution to return. This is an internal error: the Love-number "
            "method was checked to be a radial one before the solve ran.")
    # The solution takes the world with it, so its interior getters keep answering.
    world_anchor = _ProfileWorldAnchor()
    world_anchor.world_sptr = world_sptr
    solution = RadialSolverSolution._adopt(move(storage_uptr), world_anchor)
    rs_error_code = solution.error_code

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
        cy_check_surface_solve_conditioning(solution.surface_solve_amplification, c_integration_rtol)

    return solution
