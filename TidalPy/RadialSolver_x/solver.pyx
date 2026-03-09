# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Top-level radial solver entry point.
# Orchestrates EOS → shooting/matrix → love number computation.

from libc.stdio cimport printf
from libc.math cimport fabs
from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector
from libcpp.string cimport string as cpp_string

import numpy as np
cimport numpy as cnp
cnp.import_array()

from CyRK cimport PreEvalFunc, ODEMethod

from TidalPy.logger import get_logger
from TidalPy.exceptions import UnknownModelError, ArgumentException, SolutionFailedError

from TidalPy.constants cimport d_NAN, TidalPyConfig, get_shared_config_address, tidalpy_config_ptr

# Wire up the pointer at import time
if tidalpy_config_ptr == NULL:
    tidalpy_config_ptr = get_shared_config_address()

from TidalPy.utilities.math.numerics cimport c_isclose
from TidalPy.utilities.dimensions.nondimensional cimport NonDimensionalScalesCC, cf_build_nondimensional_scales

from TidalPy.RadialSolver_x.rs_solution cimport RadialSolverSolution, c_RadialSolutionStorage
from TidalPy.RadialSolver_x.shooting cimport cf_shooting_solver
from TidalPy.RadialSolver_x.matrix cimport c_matrix_propagate

# EOS Imports (Material_x)
from TidalPy.Material_x.eos cimport c_EOS_ODEInput, c_solve_eos
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.Material_x.eos.methods.interpolate cimport c_InterpolateEOSInput, c_preeval_interpolate


# EOS method constant (matches old module's EOS_INTERPOLATE_METHOD_INT)
cdef int C_EOS_INTERPOLATE_METHOD_INT = 0

log = get_logger("TidalPy")


cdef int cf_radial_solver(
        c_RadialSolutionStorage* solution_storage_ptr,
        size_t total_slices,
        double* radius_array_in_ptr,
        double* density_array_in_ptr,
        double complex* complex_bulk_modulus_in_ptr,
        double complex* complex_shear_modulus_in_ptr,
        double frequency,
        double planet_bulk_density,
        size_t num_layers,
        int* layer_types_ptr,
        int* is_static_bylayer_ptr,
        int* is_incompressible_bylayer_ptr,
        double surface_pressure,
        int degree_l,
        size_t num_bc_models,
        int* bc_models_ptr,
        int core_model,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        ODEMethod integration_method_int,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_bylayer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool nondimensionalize,
        cpp_bool use_prop_matrix,
        int* eos_integration_method_int_bylayer_ptr,
        ODEMethod eos_integration_method,
        double eos_rtol,
        double eos_atol,
        double eos_pressure_tol,
        int eos_max_iters,
        cpp_bool verbose,
        cpp_bool warnings,
        ) noexcept nogil:

    cdef size_t layer_i, slice_i

    # Figure out how many slices are in each layer
    cdef vector[size_t] first_slice_index_by_layer_vec = vector[size_t]()
    first_slice_index_by_layer_vec.resize(num_layers)

    cdef vector[size_t] num_slices_by_layer_vec = vector[size_t]()
    num_slices_by_layer_vec.resize(num_layers)

    cdef size_t layer_slices       = 0
    cdef size_t interface_check    = 0
    cdef cpp_bool top_layer        = False
    cdef double radius_check       = d_NAN
    cdef double layer_upper_radius = d_NAN

    # Pull out raw pointer
    cdef c_EOSSolution* eos_solution_storage_ptr = solution_storage_ptr.get_eos_solution_ptr()

    # Physical parameters
    cdef double radius_planet
    cdef double G_to_use                = d_NAN
    cdef double radius_planet_to_use    = d_NAN
    cdef double bulk_density_to_use     = d_NAN
    cdef double frequency_to_use        = d_NAN
    cdef double starting_radius_to_use  = d_NAN
    cdef double surface_pressure_to_use = d_NAN

    # EOS variables
    cdef size_t bottom_slice_index
    cdef vector[PreEvalFunc] eos_function_bylayer_vec = vector[PreEvalFunc]()
    eos_function_bylayer_vec.resize(num_layers)
    cdef c_EOS_ODEInput eos_input
    cdef vector[c_EOS_ODEInput] eos_inputs_bylayer_vec = vector[c_EOS_ODEInput]()
    eos_inputs_bylayer_vec.reserve(num_layers)
    cdef vector[c_InterpolateEOSInput] specific_eos_input_bylayer_vec = vector[c_InterpolateEOSInput]()
    specific_eos_input_bylayer_vec.reserve(num_layers)
    cdef char* specific_eos_char_ptr = NULL

    # Ensure there is at least one layer
    if num_layers <= 0:
        solution_storage_ptr.error_code = -5
        solution_storage_ptr.message = cpp_string(b'RadialSolver_x:: requires at least one layer, zero provided.\n')
        if verbose:
            printf(solution_storage_ptr.message.c_str())
        return solution_storage_ptr.error_code

    if solution_storage_ptr.error_code == 0:
        top_layer = False
        for layer_i in range(num_layers):
            if layer_i == num_layers - 1:
                top_layer = True

            # Determine starting slice index
            if layer_i == 0:
                first_slice_index_by_layer_vec[layer_i] = 0
            else:
                first_slice_index_by_layer_vec[layer_i] = first_slice_index_by_layer_vec[layer_i - 1] + num_slices_by_layer_vec[layer_i - 1]

            layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i]

            layer_slices = 0
            interface_check = 0
            for slice_i in range(first_slice_index_by_layer_vec[layer_i], total_slices):
                radius_check = radius_array_in_ptr[slice_i]

                if c_isclose(radius_check, layer_upper_radius, 1.0e-9, 0.0):
                    interface_check += 1
                    if interface_check > 1:
                        break
                elif radius_check > layer_upper_radius:
                    break
                layer_slices += 1

            if layer_slices < 5:
                solution_storage_ptr.error_code = -5
                solution_storage_ptr.message = cpp_string(b'RadialSolver_x:: At least five layer slices per layer are required.\n')
                if verbose:
                    printf(solution_storage_ptr.message.c_str())
                return solution_storage_ptr.error_code

            num_slices_by_layer_vec[layer_i] = layer_slices

    # Get other needed inputs
    radius_planet = radius_array_in_ptr[total_slices - 1]

    cdef NonDimensionalScalesCC non_dim_scales
    if nondimensionalize and solution_storage_ptr.error_code == 0:
        cf_build_nondimensional_scales(
            &non_dim_scales,
            frequency,
            radius_planet,
            planet_bulk_density
            )

        for slice_i in range(total_slices):
            radius_array_in_ptr[slice_i]          /= non_dim_scales.length_conversion
            density_array_in_ptr[slice_i]         /= non_dim_scales.density_conversion
            complex_bulk_modulus_in_ptr[slice_i]  /= non_dim_scales.pascal_conversion
            complex_shear_modulus_in_ptr[slice_i] /= non_dim_scales.pascal_conversion

        for layer_i in range(num_layers):
            eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i] /= non_dim_scales.length_conversion

        G_to_use                = tidalpy_config_ptr.d_G / (non_dim_scales.length3_conversion / (non_dim_scales.mass_conversion * non_dim_scales.second2_conversion))
        radius_planet_to_use    = radius_planet / non_dim_scales.length_conversion
        bulk_density_to_use     = planet_bulk_density / non_dim_scales.density_conversion
        frequency_to_use        = frequency / (1.0 / non_dim_scales.second_conversion)
        surface_pressure_to_use = surface_pressure / non_dim_scales.pascal_conversion
        starting_radius_to_use  = starting_radius / non_dim_scales.length_conversion

        solution_storage_ptr.change_radius_array(radius_array_in_ptr, total_slices, True)
    else:
        G_to_use                = tidalpy_config_ptr.d_G
        radius_planet_to_use    = radius_planet
        bulk_density_to_use     = planet_bulk_density
        frequency_to_use        = frequency
        surface_pressure_to_use = surface_pressure
        starting_radius_to_use  = starting_radius

    # Solve the equation of state
    if solution_storage_ptr.error_code == 0:
        for layer_i in range(num_layers):
            if eos_integration_method_int_bylayer_ptr[layer_i] == C_EOS_INTERPOLATE_METHOD_INT:
                eos_function_bylayer_vec[layer_i] = c_preeval_interpolate

                bottom_slice_index = first_slice_index_by_layer_vec[layer_i]

                specific_eos_input_bylayer_vec.emplace_back(
                    num_slices_by_layer_vec[layer_i],
                    &radius_array_in_ptr[bottom_slice_index],
                    &density_array_in_ptr[bottom_slice_index],
                    <double complex*>&complex_bulk_modulus_in_ptr[bottom_slice_index],
                    <double complex*>&complex_shear_modulus_in_ptr[bottom_slice_index],
                    )
                specific_eos_char_ptr = <char*>&specific_eos_input_bylayer_vec.back()

                eos_inputs_bylayer_vec.emplace_back(
                    G_to_use,
                    radius_planet_to_use,
                    specific_eos_char_ptr,
                    False,
                    False,
                    False
                )
            else:
                solution_storage_ptr.error_code = -250
                break

    if solution_storage_ptr.error_code == 0:
        c_solve_eos(
            eos_solution_storage_ptr,
            eos_function_bylayer_vec,
            eos_inputs_bylayer_vec,
            bulk_density_to_use,
            surface_pressure_to_use,
            G_to_use,
            eos_integration_method,
            eos_rtol,
            eos_atol,
            eos_pressure_tol,
            eos_max_iters,
            verbose
            )

    # Run requested radial solver method
    cdef int sub_process_error_code = 0
    if eos_solution_storage_ptr.success and solution_storage_ptr.error_code == 0:
        if use_prop_matrix:
            sub_process_error_code = c_matrix_propagate(
                solution_storage_ptr,
                frequency_to_use,
                bulk_density_to_use,
                first_slice_index_by_layer_vec.data(),
                num_slices_by_layer_vec.data(),
                num_layers,
                num_bc_models,
                bc_models_ptr,
                G_to_use,
                degree_l,
                starting_radius_to_use,
                start_radius_tolerance,
                core_model,
                verbose
                )
        else:
            sub_process_error_code = cf_shooting_solver(
                solution_storage_ptr,
                frequency_to_use,
                bulk_density_to_use,
                layer_types_ptr,
                is_static_bylayer_ptr,
                is_incompressible_bylayer_ptr,
                first_slice_index_by_layer_vec.data(),
                num_slices_by_layer_vec.data(),
                num_layers,
                num_bc_models,
                bc_models_ptr,
                G_to_use,
                degree_l,
                use_kamata,
                starting_radius_to_use,
                start_radius_tolerance,
                integration_method_int,
                integration_rtol,
                integration_atol,
                scale_rtols_bylayer_type,
                max_num_steps,
                expected_size,
                max_ram_MB,
                max_step,
                verbose
                )

    # Finalize
    if nondimensionalize:
        solution_storage_ptr.dimensionalize_data(&non_dim_scales, True)

        for slice_i in range(total_slices):
            radius_array_in_ptr[slice_i]          *= non_dim_scales.length_conversion
            density_array_in_ptr[slice_i]         *= non_dim_scales.density_conversion
            complex_bulk_modulus_in_ptr[slice_i]  *= non_dim_scales.pascal_conversion
            complex_shear_modulus_in_ptr[slice_i] *= non_dim_scales.pascal_conversion

    if solution_storage_ptr.success:
        solution_storage_ptr.find_love()

    return solution_storage_ptr.error_code


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
        cpp_bool perform_checks = True,
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

    cdef size_t layer_i, slice_i
    cdef double last_layer_r = 0.0
    cdef size_t total_slices = radius_array.size
    cdef size_t num_layers   = len(layer_types)
    cdef size_t layer_check, slice_check
    cdef cpp_bool top_layer
    cdef double last_layer_radius
    cdef double layer_radius
    cdef double radius_check, last_radius_check

    # Perform checks
    if perform_checks:
        if density_array.size != total_slices:
            raise ArgumentException("`density_array` array must be the same size as radius array.")
        if complex_bulk_modulus_array.size != total_slices:
            raise ArgumentException("`complex_bulk_modulus_array` array must be the same size as radius array.")
        if complex_shear_modulus_array.size != total_slices:
            raise ArgumentException("`complex_shear_modulus_array` array must be the same size as radius array.")

        if len(is_static_bylayer) != num_layers:
            raise ArgumentException('Number of `is_static_bylayer` must match number of `layer_types`.')
        if len(is_incompressible_bylayer) != num_layers:
            raise ArgumentException('Number of `is_incompressible_bylayer` must match number of `layer_types`.')
        if upper_radius_bylayer_array.size != num_layers:
            raise ArgumentException('Number of `upper_radius_by_layer` must match number of `layer_types`.')

        last_layer_r = 0.0
        for layer_i in range(num_layers):
            if upper_radius_bylayer_array[layer_i] <= last_layer_r:
                raise ArgumentException("`upper_radius_bylayer_array` must be in ascending order.")
            last_layer_r = upper_radius_bylayer_array[layer_i]

        if fabs(frequency) < tidalpy_config_ptr.d_MIN_FREQUENCY:
            raise ValueError('Forcing frequency is too small (are you sure you are in rad s-1?).')
        elif fabs(frequency) > tidalpy_config_ptr.d_MAX_FREQUENCY:
            raise ValueError('Forcing frequency is too large (are you sure you are in rad s-1?).')

        if use_prop_matrix:
            if num_layers > 1:
                raise NotImplementedError("Currently, TidalPy's propagation matrix technique only works for 1-layer worlds.")
            if layer_types[0].lower() != 'solid':
                raise ArgumentException("The Propagation matrix technique only works for solid layers.")
            if not is_static_bylayer[0]:
                raise ArgumentException("The Propagation matrix technique does not allow for dynamic layers.")
            if not is_incompressible_bylayer[0]:
                raise ArgumentException("The Propagation matrix technique does not allow for compressible layers.")

        if (starting_radius != 0.0) and (starting_radius > 0.90 * radius_array[total_slices - 1]):
            raise ArgumentException('Starting radius is above 90% of the planet radius. Try a lower radius.')

        if radius_array[0] != 0.0:
            raise ArgumentException('Radius array must start at zero.')

        last_layer_radius = 0.0
        for layer_i in range(num_layers):
            if layer_i == num_layers - 1:
                top_layer = True
            else:
                top_layer = False
            layer_radius = upper_radius_bylayer_array[layer_i]

            slice_check = 0
            layer_check = 0
            last_radius_check = 0.0
            for slice_i in range(total_slices):
                radius_check = radius_array[slice_i]
                if radius_check < 0.0:
                    raise ArgumentException("A negative radius value was found in `radius_array`.")
                if radius_check < last_radius_check:
                    raise ArgumentException("Radius array must be in ascending order.")
                if c_isclose(radius_check, layer_radius, 1.0e-9, 0.0):
                    layer_check += 1
                if last_layer_radius <= radius_check <= layer_radius:
                    slice_check += 1
                last_radius_check = radius_check
            last_layer_radius = layer_radius

            if slice_check < 5:
                raise ArgumentException("A minimum of 5 sub-slices (including top and bottom) are required for each layer.")

            if top_layer:
                if layer_check != 1:
                    raise ArgumentException(f"Radius of layer {layer_i} found {layer_check} times. Expected 1 time (non-interface layer).")
            else:
                if layer_check != 2:
                    raise ArgumentException(f"Radius of layer {layer_i} found {layer_check} times. Expected 2 times (interface layer).")

    # Build arrays of assumptions
    cdef vector[int] layer_types_vec = vector[int]()
    layer_types_vec.resize(num_layers)
    cdef int* layer_types_ptr = layer_types_vec.data()

    cdef vector[int] layer_assumptions_vec = vector[int]()
    layer_assumptions_vec.resize(2 * num_layers)
    cdef int* is_static_bylayer_ptr         = &layer_assumptions_vec[0]
    cdef int* is_incompressible_bylayer_ptr = &layer_assumptions_vec[num_layers]

    cdef str layer_type
    cdef cpp_bool dynamic_liquid = False

    for layer_i in range(num_layers):
        layer_type = layer_types[layer_i]

        is_static_bylayer_ptr[layer_i]         = is_static_bylayer[layer_i]
        is_incompressible_bylayer_ptr[layer_i] = is_incompressible_bylayer[layer_i]

        if not dynamic_liquid:
            if (layer_type == 1) and not is_static_bylayer_ptr[layer_i]:
                dynamic_liquid = True

        if layer_type.lower() == 'solid':
            layer_types_ptr[layer_i] = 0
        elif layer_type.lower() == 'liquid':
            layer_types_ptr[layer_i] = 1
        else:
            layer_types_ptr[layer_i] = -1
            raise UnknownModelError(f"Layer type {layer_type} is not supported. Currently supported types: 'solid', 'liquid'.")

    if perform_checks:
        if dynamic_liquid and fabs(frequency) < 2.5e-5:
            if warnings:
                log.warning(
                    'Dynamic liquid layer detected in RadialSolver for a small frequency.'
                    'Results may be unstable. Extra care is advised!'
                    )

    # Convert integration methods from string to int
    cdef str integration_method_lower = integration_method.lower()
    cdef ODEMethod integration_method_int = ODEMethod.NO_METHOD_SET
    if integration_method_lower == 'rk45':
        integration_method_int = ODEMethod.RK45
    elif integration_method_lower == 'rk23':
        integration_method_int = ODEMethod.RK23
    elif integration_method_lower == 'dop853':
        integration_method_int = ODEMethod.DOP853
    else:
        raise UnknownModelError(f"Unsupported integration method provided: {integration_method_lower}.")

    cdef str eos_integration_method_lower = eos_integration_method.lower()
    cdef ODEMethod eos_integration_method_int = ODEMethod.NO_METHOD_SET
    if eos_integration_method_lower == 'rk45':
        eos_integration_method_int = ODEMethod.RK45
    elif eos_integration_method_lower == 'rk23':
        eos_integration_method_int = ODEMethod.RK23
    elif eos_integration_method_lower == 'dop853':
        eos_integration_method_int = ODEMethod.DOP853
    else:
        raise UnknownModelError(f"Unsupported EOS integration method provided: {eos_integration_method_lower}.")

    # Convert EOS methods
    cdef str eos_method_str
    cdef vector[int] eos_integration_method_int_bylayer = vector[int]()
    eos_integration_method_int_bylayer.resize(num_layers)
    cdef int* eos_integration_method_int_bylayer_ptr = eos_integration_method_int_bylayer.data()

    if eos_method_bylayer is None:
        for layer_i in range(num_layers):
            eos_integration_method_int_bylayer[layer_i] = C_EOS_INTERPOLATE_METHOD_INT
    else:
        for layer_i in range(num_layers):
            eos_method_str = eos_method_bylayer[layer_i].lower()
            if eos_method_str == 'interpolate':
                eos_integration_method_int_bylayer[layer_i] = C_EOS_INTERPOLATE_METHOD_INT
            else:
                raise NotImplementedError("Unknown EOS method provided.")

    # Build boundary condition models
    cdef int[5] bc_models
    cdef size_t num_bc_models
    cdef int* bc_models_ptr = &bc_models[0]
    cdef str solve_for_tmp

    if solve_for is None:
        num_bc_models = 1
        bc_models_ptr[0] = 1
    else:
        if type(solve_for) != tuple:
            raise ArgumentException(
                '`solve_for` argument must be a tuple of strings. For example:\n'
                '   ("tidal",)  # If you just want tidal Love numbers.\n'
                '   ("tidal", "loading")  # If you want tidal and loading Love numbers.'
                )
        num_bc_models = len(solve_for)
        for i in range(num_bc_models):
            solve_for_tmp = solve_for[i]
            if solve_for_tmp == "free":
                bc_models_ptr[i] = 0
            elif solve_for_tmp == "tidal":
                bc_models_ptr[i] = 1
            elif solve_for_tmp == "loading":
                bc_models_ptr[i] = 2
            else:
                raise ArgumentException(f"Unsupported value provided for `solve_for`: {solve_for_tmp}.")

    # Get raw array pointers
    cdef double* radius_array_ptr  = &radius_array[0]
    cdef double* density_array_ptr = &density_array[0]
    cdef double complex* complex_shear_modulus_ptr = <double complex*>&complex_shear_modulus_array[0]
    cdef double complex* complex_bulk_modulus_ptr  = <double complex*>&complex_bulk_modulus_array[0]

    # Build solution storage
    cdef RadialSolverSolution solution = RadialSolverSolution(
        num_bc_models,
        upper_radius_bylayer_array,
        radius_array,
        degree_l
        )

    solution.set_model_names(bc_models_ptr)

    # Run radial solver
    cdef rs_error_code = 0
    rs_error_code = cf_radial_solver(
        solution.solution_storage_uptr.get(),
        total_slices,
        radius_array_ptr,
        density_array_ptr,
        complex_bulk_modulus_ptr,
        complex_shear_modulus_ptr,
        frequency,
        planet_bulk_density,
        num_layers,
        layer_types_ptr,
        is_static_bylayer_ptr,
        is_incompressible_bylayer_ptr,
        surface_pressure,
        degree_l,
        num_bc_models,
        bc_models_ptr,
        core_model,
        use_kamata,
        starting_radius,
        start_radius_tolerance,
        integration_method_int,
        integration_rtol,
        integration_atol,
        scale_rtols_bylayer_type,
        max_num_steps,
        expected_size,
        max_ram_MB,
        max_step,
        nondimensionalize,
        use_prop_matrix,
        eos_integration_method_int_bylayer_ptr,
        eos_integration_method_int,
        eos_rtol,
        eos_atol,
        eos_pressure_tol,
        eos_max_iters,
        verbose,
        warnings
        )

    # Finalize
    solution.finalize_python_storage()

    if log_info:
        solution.print_diagnostics(print_diagnostics=False, log_diagnostics=True)

    if ((not solution.success) or (rs_error_code < 0)) and raise_on_fail:
        if "not implemented" in solution.message:
            raise NotImplementedError(solution.message)
        else:
            raise SolutionFailedError(solution.message)

    if warnings:
        if np.any(solution.steps_taken > 7_000):
            log.warning(f"Large number of steps taken found in radial solver solution (max = {np.max(solution.steps_taken)}).")

    return solution
