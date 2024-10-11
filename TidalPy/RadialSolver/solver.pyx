# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport fabs, NAN
from libcpp cimport bool as cpp_bool

import numpy as np
cimport numpy as np
np.import_array()

from CyRK cimport PreEvalFunc
from CyRK.utils.vector cimport vector
from CyRK.utils.utils cimport allocate_mem, reallocate_mem, free_mem

from TidalPy.logger import get_logger
from TidalPy.exceptions import UnknownModelError
from TidalPy.utilities.dimensions.nondimensional cimport (
    cf_non_dimensionalize_physicals,
    cf_redimensionalize_physicals,
    cf_redimensionalize_radial_functions
    )
from TidalPy.RadialSolver.solutions cimport RadialSolverSolution, RadialSolutionStorageCC
from TidalPy.RadialSolver.shooting cimport cf_shooting_solver
from TidalPy.RadialSolver.matrix cimport cf_matrix_propagate

# EOS Imports (note these will change in a future release)
from TidalPy.Material.eos.interpolate cimport preeval_interpolate
from TidalPy.Material.eos.ode cimport EOS_ODEInput
from TidalPy.Material.eos.solver cimport EOSSolutionVec, solve_eos

ctypedef EOS_ODEInput* EOS_ODEInputPtr

log = get_logger(__name__)


def radial_solver(
        double[::1] radius_array,
        double[::1] density_array,
        double[::1] bulk_modulus_array,
        double complex[::1] complex_shear_modulus_array,
        double frequency,
        double planet_bulk_density,
        tuple layer_types,
        tuple is_static_by_layer,
        tuple is_incompressible_by_layer,
        tuple upper_radius_by_layer,
        double surface_pressure = 0.0,
        unsigned int degree_l = 2,
        tuple solve_for = None,
        unsigned char core_condition = 0,
        cpp_bool use_kamata = False,
        str integration_method = 'RK45',
        double integration_rtol = 1.0e-6,
        double integration_atol = 1.0e-12,
        cpp_bool scale_rtols_by_layer_type = False,
        size_t max_num_steps = 500_000,
        size_t expected_size = 500,
        size_t max_ram_MB = 500,
        double max_step = 0,
        cpp_bool limit_solution_to_radius = True,
        cpp_bool nondimensionalize = True,
        cpp_bool use_prop_matrix = False,
        cpp_bool verbose = False,
        cpp_bool warnings = True,
        cpp_bool raise_on_fail = False,
        unsigned int eos_integration_method = 2,
        double eos_rtol = 1.0e-6,
        double eos_atol = 1.0e-8,
        double eos_pressure_tol = 0.1,
        unsigned int eos_max_iters = 40
        ):
    """
    Solves the viscoelastic-gravitational problem for a planet comprised of solid and liquid layers.

    See Takeuchi and Saito (1972) for more details on this method.

    Parameters
    ----------
    radius_array : np.ndarray[dtype=np.float64]
        Radius values defined at slices throughout the planet [m].
    density_array : np.ndarray[dtype=np.float64]
        Density at each radius [kg m-3].
    bulk_modulus_array : np.ndarray[dtype=np.float64]
        Bulk modulus at each radius [Pa].
    complex_shear_modulus_array : np.ndarray[dtype=np.complex128]
        Complex shear modulus at each radius [Pa].
        This should be the result of applying a rheological model to the static shear modulus, viscosity, and
        forcing frequency.
    frequency : float64
        Forcing frequency [rad s-1]
    planet_bulk_density : float64
        Bulk density of the planet [kg m-3].
    layer_types : tuple[string, ...] (Size = number of layers)
        Indicator of layer type. Current options are:
            - "solid"
            - "liquid"
    is_static_by_layer : tuple[bool, ...] (Size = number of layers)
        Flag declaring if each layer uses the static (True) or dynamic (False) assumption.
    is_incompressible_by_layer : tuple[bool, ...] (Size = number of layers)
        Flag declaring if each layer is incompressible (True) or compressible (False).
    upper_radius_by_layer : tuple[float64, ...] (Size = number of layers)
        Tuple of the upper radius of each layer.
        Used to determine physical structure of planet.
    surface_pressure: float64, default=0
        The pressure at the surface of the planet (defined as radius_array[-1]). [Pa]
        Used for EOS calculations.
    degree_l : uint32, default=2
        Harmonic degree.
    solve_for : tuple[str, ...] (Size = number of requested solutions), default=None
        RadialSolver allows multiple solutions to be solved simultaneously, avoiding repeated integration.
        This parameter is a tuple of requested solutions. If `None`, then only "tidal" will be solved for.
        Options that are currently supported (note these are case sensitive):
            - "tidal": Tidal forcing boundary conditions.
            - "loading": Surface loading boundary conditions.
            - "free": Free surface boundary conditions.
        For example, if you want the tidal and loading solutions then you can set "solve_for=('tidal', 'loading')".
    core_condition : unsigned char, default=0
        Only used with `use_prop_matrix=True`. Defines the starting conditions at the inner boundary of the planet.
            - 0: Henning & Hurford (2014): "At the core, a special seed matrix Bcore is created with only three columns,
                 equal to the first, second, and third columns of Y for the properties at the base layer."
            - 1: Roberts & Nimmo (2008): liquid innermost zone.
            - 2: Solid innermost zone.
            - 3: Different solid innermost zone.
    use_kamata : bool, default=False
        If True, then the starting solution at the core will be based on equations from Kamata et al (2015; JGR:P)
        Otherwise, starting solution will be based on Takeuchi and Saito (1972)
    integration_method : int32, default=1
        Which CyRK integration protocol should be used. Options that are currently available are:
            - 0: Runge-Kutta 2(3)
            - 1: Runge-Kutta 4(5)
            - 2: Runge-Kutta / DOP 8(5/3)
    integration_rtol : float64, default=1.0e-4
        Relative integration tolerance. Lower tolerance will lead to more precise results at increased computation.
    integration_atol : float64, default=1.0e-12
        Absolute integration tolerance (when solution is near 0).
    scale_rtols_by_layer_type : bool, default=True
        If True, then each layer will be imparted with a different relative tolerance. Liquid layers will have a lower
        rtol which has been found to improve stability.
    max_num_steps : uint32, default=500,000
        Maximum number of integration steps allowed before solver forces a failed result.
        Useful to set to a lower number if running many calls that may be exploring an unstable parameter space for
        example, during a MCMC run.
    expected_size : uint32, default=500
        Anticipated number of integration steps that will be required per solution. Tweaking this may have a minor
        impact on performance. It is advised to leave it between 200--1000.
        Setting it too low can lead to bad performance.
    max_ram_MB : uint32, default=500
        Maximum amount of RAM in MB the integrator is allowed (ignoring some housekeeping usage).
        The integrator will use this value and the "max_num_steps" to determine a true limit on the maximum number
        steps allowed (picking the lower value.). If system RAM is limited (or if multiple RadialSolvers will run in
        parallel) it maybe worth setting this lower than the default.
        The default of 500MB is equivalent to a max_num_steps ~ 5 Mllion.
    max_step : float64, default=0
        Maximum step size the adaptive step size integrator is allowed to take. 
        Setting to 0 (default) will tell the integrator to determine an ideal max step size.
    limit_solution_to_radius : bool, default=True
        If True, then the solution will be limited to the points passed by the radius array.
    nondimensionalize : bool, default=True
        If Ture, then inputs will be non-dimensionalized before integration is performed.
        Results will be redimensionalized before being returned.
    use_prop_matrix : bool, default=False
        If True, RadialSolver will use a propagation matrix method rather than the default shooting method.
        Note that many of the parameters set by this function are only applicable to the shooting method and
        may not be passed to the propagation matrix solver.
        See more about the prop-matrix method in `TidalPy.RadialSolver.PropMatrix`.
    verbose : bool, default=False
        If True, then additioal information will be printed to the terminal during the solution. 
    warnings : bool, default=True
        If True, then warnings will be printed to the terminal during the solution. 
    raise_on_fail : bool, default=False
        If Ture, then the solver will raise an exception if integration was not successful. By default RadialSolver
        fails silently.
    eos_integration_method : unsigned int, default = 2
        Integration method used to solve for the planet's equation of state.
    eos_rtol : double, default = 1.0e-6
        Integration relative tolerance for equation of state solver.
    eos_atol : double, default = 1.0e-8
        Integration absolute tolerance for equation of state solver.
    eos_pressure_tol : double, default = 0.1
        Tolerance used when fitting to the surface pressure in the equation of state solver.
    eos_max_iters : unsigned int, default = 40
        Maximum number of iterations used to converge surface pressure in equation of state solver. 
    
    Returns
    -------
    solution : RadialSolverSolution
        Solution to the viscoelastic-gravitational problem inside a planet.
        Solution attributes:
            - solution.results : complex128, shape=(6 * num_ytypes, n_slice)
                Numerical solution throughout the planet.
            - solution.success : bool
                Flag if integration and subsequent collapse occured without error.
            - solution.message : str
                Feedback string useful for debugging.
    """
    cdef size_t i

    # Perform checks and make conversions from python to c
    cdef size_t total_slices
    total_slices = radius_array.size
    assert density_array.size               == total_slices
    assert bulk_modulus_array.size          == total_slices
    assert complex_shear_modulus_array.size == total_slices

    # Unpack inefficient user-provided tuples into bool arrays and pass by pointer
    cdef size_t num_layers
    num_layers = len(layer_types)

    # Check that number of assumptions match.
    if len(is_static_by_layer) != num_layers:
        raise AttributeError('Number of `is_static_by_layer` must match number of `layer_types`.')
    if len(is_incompressible_by_layer) != num_layers:
        raise AttributeError('Number of `is_incompressible_by_layer` must match number of `layer_types`.')
    if len(upper_radius_by_layer) != num_layers:
        raise AttributeError('Number of `upper_radius_by_layer` must match number of `layer_types`.')

    # Build array of assumptions
    # OPT: Perhaps set a maximum number of layers then we can put these on the stack rather than heap allocating them.
    cdef int* layer_assumptions_ptr = <int *> allocate_mem(
        3 * num_layers * sizeof(int),
        'layer_assumptions_ptr (radial_solver; init)'
        )
    cdef int* layer_types_ptr                = &layer_assumptions_ptr[0]
    cdef int* is_static_by_layer_ptr         = &layer_assumptions_ptr[num_layers]
    cdef int* is_incompressible_by_layer_ptr = &layer_assumptions_ptr[2 * num_layers]
    cdef double* upper_radius_by_layer_ptr = <double *> allocate_mem(
        num_layers * sizeof(double),
        'upper_radius_by_layer_ptr (radial_solver; init)'
        )
    
    cdef str layer_type
    cdef bint dynamic_liquid = False

    # Pull out information for each layer and store in heap memory
    for i in range(num_layers):
        layer_type                        = layer_types[i]
        is_static_by_layer_ptr[i]         = is_static_by_layer[i]
        is_incompressible_by_layer_ptr[i] = is_incompressible_by_layer[i]
        upper_radius_by_layer_ptr[i]      = upper_radius_by_layer[i]

        if not dynamic_liquid:
            if (layer_type == 1) and not is_static_by_layer_ptr[i]:
                # There is at least one dynamic liquid layer
                dynamic_liquid = True

        # Convert user-provided strings to ints for the layer type
        if layer_type.lower() == 'solid':
            layer_types_ptr[i] = 0
        elif layer_type.lower() == 'liquid':
            layer_types_ptr[i] = 1
        else:
            layer_types_ptr[i] = -1
            log.error(f"Layer type {layer_type} is not supported. Currently supported types: 'solid', 'liquid'.")
            raise UnknownModelError(f"Layer type {layer_type} is not supported. Currently supported types: 'solid', 'liquid'.")
    
    # Check for dynamic liquid layer stability
    if dynamic_liquid and fabs(frequency) < 2.5e-5:
        # TODO: check that this frequency is a decent cutoff (based on a 3 day period).
        #    Initial work suggests that low density layers do not suffer from the same instability problems.
        # TODO: Add density or combo factor to better indicate when a solution is likely to be unstable?
        # See Issue #55
        if warnings:
            log.warning(
                'Dynamic liquid layer detected in RadialSolver for a small frequency.'
                'Results may be unstable. Extra care is advised!'
                )
    
    # Convert integration method to int
    cdef str integration_method_lower = integration_method.lower()
    cdef unsigned char integration_method_int
    if integration_method_lower == 'rk45':
        integration_method_int = 1
    elif integration_method_lower == 'rk23':
        integration_method_int = 0
    elif integration_method_lower == 'dop853':
        integration_method_int = 2
    else:
        log.error(f"Unsupported integration method provided: {integration_method_lower}.")
        raise UnknownModelError(f"Unsupported integration method provided: {integration_method_lower}.")
    
    # Clean up what values the solver is solving for.
    cdef int[5] bc_models
    cdef size_t num_bc_models
    cdef int* bc_models_ptr = &bc_models[0]
    cdef str solve_for_tmp 

    if solve_for is None:
        # `solve_for` was not provided. Assume the user wants the tidal, and only the tidal, solution.
        num_bc_models = 1
        bc_models_ptr[0] = 1
    else:
        if type(solve_for) != tuple:
            raise AttributeError(
                '`solve_for` argument must be a tuple of strings. For example:\n'
                '   ("tidal",)  # If you just want tidal Love numbers.\n'
                '   ("tidal", "loading")  # If you just want tidal and loading Love numbers.'
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
                raise AttributeError(f"Unsupported value provided for `solve_for`: {solve_for_tmp}.")
    
    # Get other inputs needed
    cdef double radius_planet = radius_array[total_slices - 1]

    # Build solution storage
    cdef RadialSolverSolution solution = RadialSolverSolution(total_slices, num_bc_models)
    solution.set_model_names(bc_models_ptr)
    cdef RadialSolutionStorageCC* solution_storage_ptr = solution.solution_storage_ptr

    # Convert complex-valued arrays to C++ complex pointers
    cdef double complex* complex_shear_modulus_ptr = <double complex*>&complex_shear_modulus_array[0]

    # Solve the equaiton of state for the planet

    # Non-dimensionalize inputs
    cdef cpp_bool already_nondimed = False
    cdef double G_to_use = NAN
    cdef double radius_planet_to_use = NAN
    cdef double bulk_density_to_use = NAN
    cdef double frequency_to_use = NAN
    if nondimensionalize:
        cf_non_dimensionalize_physicals(
            total_slices,
            frequency,
            radius_planet,
            planet_bulk_density,
            &radius_array[0],
            &density_array[0],
            &bulk_modulus_array[0],
            complex_shear_modulus_ptr,
            &radius_planet_to_use,
            &bulk_density_to_use,
            &frequency_to_use,
            &G_to_use
            )
        already_nondimed = True

    # TODO: For now there is only one accepted EOS, the interpolated kind. In the future additional EOS will be supplied
    # either via arguments to this function or a more OOP approach where they are built into the layers.
    # Build arrays of EOS inputs.
    cdef vector[preeval_interpolate] eos_function_bylayer_vec = vector[preeval_interpolate](0)
    eos_function_bylayer_vec.reserve(num_layers)

    # Build vector of inputs
    cdef vector[EOS_ODEInput] eos_inputs_bylayer_vec = vector[EOS_ODEInput](0)
    cdef vector[EOS_ODEInputPtr] eos_inputs_ptrs_bylayer_vec = vector[EOS_ODEInputPtr](0)
    eos_inputs_ptrs_bylayer_vec.reserve(num_layers)
    eos_inputs_bylayer_vec.reserve(num_layers)
    for i in range(num_layers):
        # TODO: For now we are only storing the interpolate version of the EOS for each layer.
        eos_function_bylayer_vec.push_back(preeval_interpolate)
        eos_inputs_ptrs_bylayer_vec.push_back(&eos_function_bylayer_vec[i])

    # Make pointers to pre-eval data
    cdef PreEvalFunc* eos_function_bylayer_ptrs = &eos_function_bylayer_vec[0]
    cdef EOS_ODEInput** eos_input_bylayer_ptrs  = &eos_inputs_ptrs_bylayer_vec[0]

    cdef EOSSolutionVec eos_result = solve_eos(
        &radius_array[0],
        total_slices,
        upper_radius_by_layer_ptr,
        num_layers,
        eos_function_bylayer_ptrs,
        eos_input_bylayer_ptrs,
        planet_bulk_density,
        surface_pressure,
        G_to_use,
        eos_integration_method,
        eos_rtol,
        eos_atol,
        eos_pressure_tol,
        eos_max_iters
        )

    # Run requested radial solver method
    try:
        if use_prop_matrix:
            cf_matrix_propagate(
                solution_storage_ptr,
                total_slices,
                &radius_array[0],
                &density_array[0],
                &bulk_modulus_array[0],
                complex_shear_modulus_ptr,
                frequency,
                planet_bulk_density,
                # TODO: In the future the propagation matrix should take in layer types and multiple layers
                # size_t num_layers,
                # int* layer_types_ptr,
                # int* is_static_by_layer_ptr,
                # int* is_incompressible_by_layer_ptr,
                # double* upper_radius_by_layer_ptr,
                num_bc_models,
                bc_models_ptr,
                degree_l,
                core_condition,
                nondimensionalize,
                verbose,
                raise_on_fail,
                already_nondimed)
        else:
            cf_shooting_solver(
                solution_storage_ptr,
                total_slices,
                &radius_array[0],
                &density_array[0],
                &gravity_array[0],
                &bulk_modulus_array[0],
                complex_shear_modulus_ptr,
                frequency,
                planet_bulk_density,
                num_layers,
                layer_types_ptr,
                is_static_by_layer_ptr,
                is_incompressible_by_layer_ptr,
                upper_radius_by_layer_ptr,
                num_bc_models,
                bc_models_ptr,
                degree_l,
                use_kamata,
                integration_method_int,
                integration_rtol,
                integration_atol,
                scale_rtols_by_layer_type,
                max_num_steps,
                expected_size,
                max_ram_MB,
                max_step,
                limit_solution_to_radius,
                nondimensionalize,
                verbose,
                raise_on_fail,
                already_nondimed)
    finally:
        # Release heap memory
        if not (layer_assumptions_ptr is NULL):
            layer_types_ptr                = NULL
            is_static_by_layer_ptr         = NULL
            is_incompressible_by_layer_ptr = NULL
            free_mem(layer_assumptions_ptr)
            layer_assumptions_ptr = NULL
        if not (upper_radius_by_layer_ptr is NULL):
            free_mem(upper_radius_by_layer_ptr)
            upper_radius_by_layer_ptr = NULL
    
    return solution
