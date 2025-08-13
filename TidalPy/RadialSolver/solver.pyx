# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False


from libc.stdio cimport printf
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.math cimport fabs
from libc.string cimport strcpy
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

from CyRK cimport PreEvalFunc

from TidalPy.logger import get_logger
from TidalPy.exceptions import UnknownModelError, ArgumentException, SolutionFailedError

from TidalPy.constants cimport d_G, d_MIN_FREQUENCY, d_MAX_FREQUENCY, d_NAN_DBL
from TidalPy.utilities.math.numerics cimport cf_isclose
from TidalPy.utilities.dimensions.nondimensional cimport NonDimensionalScalesCC, cf_build_nondimensional_scales
from TidalPy.RadialSolver.rs_solution cimport RadialSolverSolution
from TidalPy.RadialSolver.shooting cimport cf_shooting_solver
from TidalPy.RadialSolver.matrix cimport cf_matrix_propagate

# EOS Imports (note these will change in a future release)
from TidalPy.Material.eos cimport EOS_ODEInput, solve_eos
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.Material.eos.methods cimport InterpolateEOSInput, preeval_interpolate, EOS_INTERPOLATE_METHOD_INT


ctypedef EOS_ODEInput* EOS_ODEInputPtr

log = get_logger("TidalPy")


cdef int cf_radial_solver(
        shared_ptr[RadialSolutionStorageCC] solution_storage_sptr,
        size_t total_slices,
        double* radius_array_in_ptr,
        double* density_array_in_ptr,
        double complex* complex_bulk_modulus_in_ptr,
        double complex* complex_shear_modulus_in_ptr,
        double frequency,
        double planet_bulk_density,
        size_t num_layers,
        int* layer_types_ptr,
        bint* is_static_bylayer_ptr,
        bint* is_incompressible_bylayer_ptr,
        double surface_pressure,
        int degree_l,
        size_t num_bc_models,
        int* bc_models_ptr,
        int core_model,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        int integration_method_int,
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
        int eos_integration_method,
        double eos_rtol,
        double eos_atol,
        double eos_pressure_tol,
        int eos_max_iters,
        cpp_bool verbose,
        cpp_bool warnings,
        ) noexcept nogil:

    cdef size_t layer_i, slice_i
    # Feedback
    cdef char[256] message
    cdef char* message_ptr = &message[0]
    
    # Figure out how many slices are in each layer
    cdef vector[size_t] first_slice_index_by_layer_vec = vector[size_t]()
    first_slice_index_by_layer_vec.resize(num_layers)

    cdef vector[size_t] num_slices_by_layer_vec = vector[size_t]()
    num_slices_by_layer_vec.resize(num_layers)

    cdef size_t layer_slices       = 0
    cdef size_t interface_check    = 0
    cdef cpp_bool top_layer        = False
    cdef double radius_check       = d_NAN_DBL
    cdef double layer_upper_radius = d_NAN_DBL

    # Pull out raw pointers to avoid repeated calls to the getter
    cdef RadialSolutionStorageCC* solution_storage_ptr = solution_storage_sptr.get()
    cdef EOSSolutionCC* eos_solution_storage_ptr       = solution_storage_ptr.get_eos_solution_ptr()

    # Physical parameters
    cdef double radius_planet
    cdef double G_to_use                = d_NAN_DBL
    cdef double radius_planet_to_use    = d_NAN_DBL
    cdef double bulk_density_to_use     = d_NAN_DBL
    cdef double frequency_to_use        = d_NAN_DBL
    cdef double starting_radius_to_use  = d_NAN_DBL
    cdef double surface_pressure_to_use = d_NAN_DBL
    
    # Equation of state variables
    cdef size_t bottom_slice_index
    cdef vector[PreEvalFunc] eos_function_bylayer_vec = vector[PreEvalFunc]()
    eos_function_bylayer_vec.resize(num_layers)
    # Build vector of EOS inputs
    cdef EOS_ODEInput eos_input
    cdef vector[EOS_ODEInput] eos_inputs_bylayer_vec = vector[EOS_ODEInput]()
    eos_inputs_bylayer_vec.reserve(num_layers)
    # EOS Input is stored as void pointers because each model may have different input requirements.
    cdef vector[InterpolateEOSInput] specific_eos_input_bylayer_vec = vector[InterpolateEOSInput]() # TODO: Only for specific EOS model; have a vector[void] ??
    specific_eos_input_bylayer_vec.reserve(num_layers)
    cdef char* specific_eos_char_ptr = NULL

    # Ensure there is at least one layer.
    if num_layers <= 0:
        solution_storage_ptr.error_code = -5
        strcpy(message_ptr, 'RadialSolver:: requires at least one layer, zero provided.\n')
        solution_storage_ptr.set_message(message_ptr)
        if verbose:
            printf(message_ptr)
        return solution_storage_ptr.error_code

    if solution_storage_ptr.error_code == 0:
        top_layer = False
        for layer_i in range(num_layers):
            if layer_i == num_layers - 1:
                top_layer = True

            # Determine starting slice index in this layer
            if layer_i == 0:
                first_slice_index_by_layer_vec[layer_i] = 0
            else:
                first_slice_index_by_layer_vec[layer_i] = first_slice_index_by_layer_vec[layer_i - 1] + num_slices_by_layer_vec[layer_i - 1]
            
            # Determine number of slices in layer
            layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i]

            layer_slices = 0
            interface_check = 0
            for slice_i in range(first_slice_index_by_layer_vec[layer_i], total_slices):
                radius_check = radius_array_in_ptr[slice_i]
                
                # TidalPy requires that each layer's upper radius be provided twice for interface layers and once for top-most layers. 
                if cf_isclose(radius_check, layer_upper_radius):
                    # Found slice that matches this layer's upper radius. We want to grab one for sure.
                    interface_check += 1
                    # We do not want to grab a second (there would be two for interface layers)
                    if interface_check > 1:
                        break
                elif radius_check > layer_upper_radius:
                    # We have passed this layer.
                    break
                layer_slices += 1
            
            if layer_slices < 5:
                solution_storage_ptr.error_code == -5
                strcpy(message_ptr, 'RadialSolver:: At least five layer slices per layer are required. Try using more slices in the input arrays.\n')
                solution_storage_ptr.set_message(message_ptr)
                if verbose:
                    printf(message_ptr)
                return solution_storage_ptr.error_code

            num_slices_by_layer_vec[layer_i] = layer_slices

    # Get other needed inputs
    radius_planet = radius_array_in_ptr[total_slices - 1]

    cdef NonDimensionalScalesCC non_dim_scales
    if nondimensionalize and solution_storage_ptr.error_code == 0:

        # Create scales used to non-dimensionalize various properties
        cf_build_nondimensional_scales(
            &non_dim_scales,
            frequency,
            radius_planet,
            planet_bulk_density
            )

        # Convert array pointers
        for slice_i in range(total_slices):
            radius_array_in_ptr[slice_i]          /= non_dim_scales.length_conversion
            density_array_in_ptr[slice_i]         /= non_dim_scales.density_conversion
            complex_bulk_modulus_in_ptr[slice_i]  /= non_dim_scales.pascal_conversion
            complex_shear_modulus_in_ptr[slice_i] /= non_dim_scales.pascal_conversion
        
        for layer_i in range(num_layers):
            eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i] /= non_dim_scales.length_conversion
            
        # Convert non-array constants
        G_to_use                = d_G / (non_dim_scales.length3_conversion / (non_dim_scales.mass_conversion * non_dim_scales.second2_conversion))
        radius_planet_to_use    = radius_planet / non_dim_scales.length_conversion
        bulk_density_to_use     = planet_bulk_density / non_dim_scales.density_conversion
        frequency_to_use        = frequency / (1. / non_dim_scales.second_conversion)
        surface_pressure_to_use = surface_pressure / non_dim_scales.pascal_conversion
        starting_radius_to_use  = starting_radius / non_dim_scales.length_conversion

        # Scale tolerances
        # TODO: What about rtol and atol for EOS and radial solver?
        # TODO How about these other tolerances?
        # eos_pressure_tol /= non_dim_scales.pascal_conversion
        # start_radius_tolerance /= non_dim_scales.length_conversion

        # Update the radius array inside the C++ classes
        solution_storage_ptr.change_radius_array(radius_array_in_ptr, total_slices, True)

    else:
        G_to_use                = d_G
        radius_planet_to_use    = radius_planet
        bulk_density_to_use     = planet_bulk_density
        frequency_to_use        = frequency
        surface_pressure_to_use = surface_pressure
        starting_radius_to_use  = starting_radius

    # Solve the equation of state for the planet
    # TODO: For now there is only one accepted EOS, the interpolated kind. In the future additional EOS will be supplied
    # either via arguments to this function or a more OOP approach where they are built into the layers.
    # Build arrays of EOS inputs.
    # Build Equation of State functions and input data structures. Record memory addresses for use by the EOS solver
    if solution_storage_ptr.error_code == 0:
        for layer_i in range(num_layers):
            # TODO: Below is specific to interpolate EOS. For now we are only storing the interpolate version of the EOS for each layer.

            if eos_integration_method_int_bylayer_ptr[layer_i] == EOS_INTERPOLATE_METHOD_INT:

                eos_function_bylayer_vec[layer_i] = preeval_interpolate

                # Build EOS input
                bottom_slice_index = first_slice_index_by_layer_vec[layer_i]

                # Build EOS input for specific EOS model
                # TODO: Below is specific to interpolate EOS.
                specific_eos_input_bylayer_vec.emplace_back(
                    num_slices_by_layer_vec[layer_i],                   # Number of slices for this layer [size_t]
                    &radius_array_in_ptr[bottom_slice_index],           # Radius array pointer [double*]
                    &density_array_in_ptr[bottom_slice_index],          # Density array pointer [double*]
                    &complex_bulk_modulus_in_ptr[bottom_slice_index],   # Complex bulk array pointer [double complex*]
                    &complex_shear_modulus_in_ptr[bottom_slice_index],  # Complex shear array pointer [double complex*]
                    )
                specific_eos_char_ptr = <char*>&specific_eos_input_bylayer_vec.back()
                
                # Build input for generalized EOS solver
                eos_inputs_bylayer_vec.emplace_back(
                    G_to_use,                 # Gravitational constant [double]
                    radius_planet_to_use,     # Planet radius [double]
                    specific_eos_char_ptr,    # void-casted pointer to model specific input (created just above)
                    False,                    # Final solve flag [bool] (will be updated by EOS solver)
                    False,                    # Final update shear flag [bool] (will be updated by EOS solver)
                    False                     # Final update bulk flag [bool] (will be updated by EOS solver)
                )
                # TODO: update bulk/shear flags are overwritten by EOS solver regardless of layer type. They should not ever need to be updated, for example shear for liquid layers?
            else:
                # Unknown EOS method
                solution_storage_ptr.error_code = -250
                break

    if solution_storage_ptr.error_code == 0:
        solve_eos(
            eos_solution_storage_ptr,   # Equation of state storage C++ class
            eos_function_bylayer_vec,   # EOS specific model function by layer array pointer [vector<PreEvalFunc>]
            eos_inputs_bylayer_vec,     # Pointer to array of EOS input pointers for each layer [vector<shared_ptr<EOS_ODEInput>>]
            bulk_density_to_use,        # Planet bulk density [double]
            surface_pressure_to_use,    # Planet surface pressure [double]
            G_to_use,                   # Gravitational constant [double]
            eos_integration_method,     # Integration method [unsigned int]
            eos_rtol,                   # Integration relative tolerance [double]
            eos_atol,                   # Integration absolute tolerance [double]
            eos_pressure_tol,           # Pressure tolerance (used for EOS convergence with surface pressure) [double]
            eos_max_iters,              # Maximum iterations to find convergence [unsigned int]
            verbose                     # Verbose flag [cpp_bool]
            )

    # Step through the radial steps to find EOS-dependent parameters
    cdef size_t i
    cdef int sub_process_error_code = 0
    if eos_solution_storage_ptr.success and solution_storage_ptr.error_code == 0:
        # Run requested radial solver method
        if use_prop_matrix:
            sub_process_error_code = cf_matrix_propagate(
                solution_storage_ptr,           # (Modified) Final radial solution storage struct pointer [RadialSolutionStorageCC*]
                frequency_to_use,               # Forcing frequency [double]
                bulk_density_to_use,            # Planet bulk density [double]
                # TODO: In the future the propagation matrix should take in layer types and multiple layers
                # int* layer_types_ptr,
                # int* is_static_bylayer_ptr,
                # int* is_incompressible_bylayer_ptr,
                first_slice_index_by_layer_vec, # First radial slice of each layer array pointer [size_t*]
                num_slices_by_layer_vec,        # Number of radial slices in each layer array pointer [size_t*]
                num_bc_models,                  # Number of boundary conditions requested by user [size_t]
                bc_models_ptr,                  # Boundary condition model int array pointer [int*]
                G_to_use,                       # Gravitational constant [double]
                degree_l,                       # Harmonic degree [unsigned int]
                starting_radius_to_use,         # Starting radius for solver. For higher degree solutions you generally want to start higher up in the planet. [double]
                start_radius_tolerance,         # Tolerance used if `starting_radius` is not provided. [double]                
                core_model,                     # Starting condition model int at the inner boundary (usually a core) see TidalPy.RadialSolver.matrix.pyx for options [unsigned char]
                verbose,                        # Verbose flag [cpp_bool]
                )
        else:
            sub_process_error_code = cf_shooting_solver(
                solution_storage_ptr,           # (Modified) Final radial solution storage struct pointer [RadialSolutionStorageCC*]
                frequency_to_use,               # Forcing frequency [double]
                bulk_density_to_use,            # Planet bulk density [double]
                layer_types_ptr,                # Layer type int  array pointer [int*]
                is_static_bylayer_ptr,          # Layer is_static flag array pointer [int*]
                is_incompressible_bylayer_ptr,  # Pointer array of layer is_incompressible flag array pointer [int*]
                first_slice_index_by_layer_vec, # First radial slice of each layer array pointer [size_t*]
                num_slices_by_layer_vec,        # Number of radial slices in each layer array pointer [size_t*]
                num_bc_models,                  # Number of boundary conditions requested by user [size_t]
                bc_models_ptr,                  # Boundary condition model int array pointer [int*]
                G_to_use,                       # Gravitational constant [double]
                degree_l,                       # Harmonic degree [unsigned int]
                use_kamata,                     # Flag to use Kamata+ (2015)'s starting conditions vs. Takeuchi+Saito (1972) [cpp_bool]
                starting_radius_to_use,         # Starting radius for solver. For higher degree solutions you generally want to start higher up in the planet. [double]
                start_radius_tolerance,         # Tolerance used if `starting_radius` is not provided. [double]
                integration_method_int,         # Integration method int (0=RK23, 1=RK45, 2=DOP853) [unsigned char]
                integration_rtol,               # Integration relative tolerance [double]
                integration_atol,               # Integration absolute tolerance [double]
                scale_rtols_bylayer_type,       # Flag for if tolerances should vary with layer type (using pre-defined scaling) [cpp_bool]
                max_num_steps,                  # Maximum number of integration steps allowed [size_t]
                expected_size,                  # Expected number of integration steps required for the average layer [size_t]
                max_ram_MB,                     # Maximum amount of ram allowed for each layer's integration (note if parallelized then radial solver will exceed this value; there is also overhead of other functions) [size_t]
                max_step,                       # Maximum allowed step size per layer [double]
                verbose,                        # Verbose flag [cpp_bool]
                )

    # Finalize solution storage
    cdef size_t solver_i
    cdef size_t bc_stride
    cdef size_t solver_stride
    cdef double complex* full_solution_ptr = NULL
    if nondimensionalize:
        # Redimensionalize user-provided inputs that are provided as pointers so that this function returns to the same state.
        
        solution_storage_sptr.get().dimensionalize_data(&non_dim_scales, True)

        # Redimensionalize user-provided inputs that are provided as pointers so that this function returns to the same state.
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

    See Takeuchi and Saito (1972) for more details on this method.

    Parameters
    ----------
    radius_array : np.ndarray[dtype=np.float64]
        Radius values defined at slices throughout the planet [m].
        It must start at zero and end at the planet's surface.
    density_array : np.ndarray[dtype=np.float64]
        Density at each radius [kg m-3].
    complex_bulk_modulus_array : np.ndarray[dtype=np.complex128]
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
    is_static_bylayer : tuple[bool, ...] (Size = number of layers)
        Flag declaring if each layer uses the static (True) or dynamic (False) assumption.
    is_incompressible_bylayer : tuple[bool, ...] (Size = number of layers)
        Flag declaring if each layer is incompressible (True) or compressible (False).
    upper_radius_by_layer : tuple[float64, ...] (Size = number of layers)
        Tuple of the upper radius of each layer.
        Used to determine physical structure of planet.
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
    starting_radius : float64, default=0.0
        The initial radius at which to start solving the viscoelastic-gravitational equations.
        For stability purposes, the higher your degree l, the higher you want your starting r. 
        If set to 0.0, the default, TidalPy will determine a good starting l based on the planet radius and degree l.
        If you provide a non-zero value then it must be in [m]
    start_radius_tolerance : float64, default=1.0e-5
        If `starting_radius` is not provided then TidalPy will use the formula:
            r = planet_radius * start_radius_tolerance^(1.0 / degree_l)
        Depending on your layer size and the degree_l, this may be in a layer above the innermost core. 
        TidalPy will perform a full planet equation of state calculation, but will skip the viscoelastic calculations
        for layers and radii below the starting value.
    nondimensionalize : bool, default=True
        If True, then inputs will be non-dimensionalized before integration (either shooting method or prop matrix) is performed.
        Results will be redimensionalized before being returned.
    use_kamata : bool, default=False
        If True, then the starting solution at the core will be based on equations from Kamata et al (2015; JGR:P)
        Otherwise, starting solution will be based on Takeuchi and Saito (1972)
    integration_method : int32, default="DOP853"
        Which CyRK integration protocol should be used. Options that are currently available are:
            - 0: Runge-Kutta 2(3)
            - 1: Runge-Kutta 4(5)
            - 2: Runge-Kutta / DOP 8(5/3)
    integration_rtol : float64, default=1.0e-5
        Relative integration tolerance. Lower tolerance will lead to more precise results at increased computation.
    integration_atol : float64, default=1.0e-8
        Absolute integration tolerance (when solution is near 0).
    scale_rtols_bylayer_type : bool, default=True
        If True, then each layer will be imparted with a different relative tolerance. Liquid layers will have a lower
        rtol which has been found to improve stability.
    max_num_steps : uint32, default=500,000
        Maximum number of integration steps allowed before solver forces a failed result.
        Useful to set to a lower number if running many calls that may be exploring an unstable parameter space for
        example, during a MCMC run.
    expected_size : uint32, default=1000
        Anticipated number of integration steps that will be required per solution. Tweaking this may have a minor
        impact on performance. It is advised to leave it between 200--1000.
        Setting it too low can lead to bad performance.
    max_ram_MB : uint32, default=500
        Maximum amount of RAM in MB the integrator is allowed (ignoring some housekeeping usage).
        The integrator will use this value and the "max_num_steps" to determine a true limit on the maximum number
        steps allowed (picking the lower value.). If system RAM is limited (or if multiple RadialSolvers will run in
        parallel) it maybe worth setting this lower than the default.
        The default of 500MB is equivalent to a max_num_steps ~ 5 Million.
    max_step : float64, default=0
        Maximum step size the adaptive step size integrator is allowed to take. 
        Setting to 0 (default) will tell the integrator to determine an ideal max step size.
    use_prop_matrix : bool, default=False
        If True, RadialSolver will use a propagation matrix method rather than the default shooting method.
        Note that many of the parameters set by this function are only applicable to the shooting method and
        may not be passed to the propagation matrix solver.
        See more about the prop-matrix method in `TidalPy.RadialSolver.PropMatrix`.
    core_model : uint32, default=0
        Only used with `use_prop_matrix=True`. Defines the starting conditions at the inner boundary of the planet.
            - 0: Henning & Hurford (2014): "At the core, a special seed matrix Bcore is created with only three columns,
                 equal to the first, second, and third columns of Y for the properties at the base layer."
            - 1: Roberts & Nimmo (2008): liquid innermost zone.
            - 2: Solid innermost zone.
            - 3: Different solid innermost zone.
            - 4: Sabadini & Veermeerson (2004), More complex interface matrix
    eos_method_by_layer : tuple, default = None
        Tuple of EOS methods for each layer. This is a tuple of strings.
        If `None` then will use default for each layer (interpolation)
        Currently supported methods:
            - "interpolation"
    surface_pressure: float64, default=0
        The pressure at the surface of the planet (defined as radius_array[-1]). [Pa]
        Used for EOS calculations.
    eos_integration_method : unsigned int, default = "DOP853"
        Integration method used to solve for the planet's equation of state.
    eos_rtol : double, default = 1.0e-3
        Integration relative tolerance for equation of state solver.
    eos_atol : double, default = 1.0e-5
        Integration absolute tolerance for equation of state solver.
    eos_pressure_tol : double, default = 0.1
        Tolerance used when fitting to the surface pressure in the equation of state solver.
    eos_max_iters : int, default = 40
        Maximum number of iterations used to converge surface pressure in equation of state solver.
    verbose : bool, default=False
        If True, then additional information will be printed to the terminal during the solution. 
    warnings : bool, default=True
        If True, then warnings will be printed to the terminal during the solution. 
    raise_on_fail : bool, default=False
        If True, then the solver will raise an exception if integration was not successful. By default RadialSolver
        fails silently.
    perform_checks : bool, default=True
        Performs sanity checks that raise python exceptions. If turned off then these checks will be skipped providing 
        some boost to performance but at the risk of uncaught exceptions (crashes).
    log_info : bool, default=False
        Flag to turn on logging of key information (diagnostic and physical) from the RadialSolverSolution.
        Note there is a performance hit if this is true, particularly if file logging is enabled.
    
    Returns
    -------
    solution : RadialSolverSolution
        Solution to the viscoelastic-gravitational problem inside a planet.
        Also contains the EOS solver solution for the entire planet.
    """

    cdef size_t layer_i, slice_i
    cdef double last_layer_r = 0.
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

        # Check that number of assumptions match.
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

        if fabs(frequency) < d_MIN_FREQUENCY:
            raise ValueError('Forcing frequency is too small (are you sure you are in rad s-1?).')
        elif fabs(frequency) > d_MAX_FREQUENCY:
            raise ValueError('Forcing frequency is too large (are you sure you are in rad s-1?).')
        
        if use_prop_matrix:
            if num_layers > 1:
                raise NotImplementedError("Currently, TidalPy's propagation matrix technique only works for 1-layer worlds. For 2 layer worlds where the lower layer is a liquid: you can start the solver at the bottom of the upper solid layer.")
            if layer_types[0].lower() != 'solid':
                raise ArgumentException("The Propagation matrix technique only works for solid layers. For liquid layers you can set layer type to solid and use a small shear modulus to mimic liquid layers.")
            if not is_static_bylayer[0]:
                raise ArgumentException("The Propagation matrix technique does not allow for dynamic layers.")
            if not is_incompressible_bylayer[0]:
                raise ArgumentException("The Propagation matrix technique does not allow for compressible layers.")
        
        if (starting_radius != 0.0) and (starting_radius > 0.90 * radius_array[total_slices - 1]):
            raise ArgumentException('Starting radius is above 90\% of the planet radius. Try a lower radius.')

        # Check radius array. TidalPy requires very specific requirements for the radius array format.
        #  1) Must start at 0.
        #  2) Radius must be ordered in ascending order.
        #  3) Each layer must have the starting radius and the ending radius. Yes that means there will be duplicate values of r at interfaces.
        #  4) There must be at least 5 sub slices in each layer
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
                    # Array must be ascending order. Duplicates are required at interfaces but nothing should be less than previous.
                    raise ArgumentException("Radius array must be in ascending order.")

                if cf_isclose(radius_check, layer_radius):
                    layer_check += 1
                if last_layer_radius <= radius_check <= layer_radius:
                    # Inside the layer.
                    slice_check += 1
                last_radius_check = radius_check
            last_layer_radius = layer_radius
            
            if slice_check < 5:
                raise ArgumentException("A minimum of 5 sub-slices (including top and bottom) are required for each layer.")

            if top_layer:
                # Top layer works for both single layer planets or multi-layer since a single layer is the top layer.
                if layer_check != 1:
                    raise ArgumentException(f"Radius of layer {layer_i} found {layer_check} times. Expected 1 time (non-interface layer).")
            else:
                if layer_check != 2:
                    raise ArgumentException(f"Radius of layer {layer_i} found {layer_check} times. Expected 2 times (interface layer).")
        
    # Build array of assumptions
    # OPT: Perhaps set a maximum number of layers then we can put these on the stack rather than heap allocating them.
    cdef vector[int] layer_types_vec = vector[int]()
    layer_types_vec.resize(num_layers)
    cdef int* layer_types_ptr = layer_types_vec.data()

    cdef vector[bint] layer_assumptions_vec = vector[bint]()
    layer_assumptions_vec.resize(2 * num_layers)
    cdef bint* is_static_bylayer_ptr         = &layer_assumptions_vec[0]
    cdef bint* is_incompressible_bylayer_ptr = &layer_assumptions_vec[num_layers]

    cdef str layer_type
    cdef cpp_bool dynamic_liquid = False

    # Pull out information for each layer and store in heap memory
    for layer_i in range(num_layers):
        layer_type = layer_types[layer_i]
        
        is_static_bylayer_ptr[layer_i]         = is_static_bylayer[layer_i]
        is_incompressible_bylayer_ptr[layer_i] = is_incompressible_bylayer[layer_i]

        if not dynamic_liquid:
            if (layer_type == 1) and not is_static_bylayer_ptr[layer_i]:
                # There is at least one dynamic liquid layer
                dynamic_liquid = True

        # Convert user-provided strings to ints for the layer type
        if layer_type.lower() == 'solid':
            layer_types_ptr[layer_i] = 0
        elif layer_type.lower() == 'liquid':
            layer_types_ptr[layer_i] = 1
        else:
            layer_types_ptr[layer_i] = -1
            raise UnknownModelError(f"Layer type {layer_type} is not supported. Currently supported types: 'solid', 'liquid'.")
    
    # Check for dynamic liquid layer stability
    if perform_checks:
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
    
    # Convert integration methods from string to int
    cdef str integration_method_lower = integration_method.lower()
    cdef int integration_method_int
    if integration_method_lower == 'rk45':
        integration_method_int = 1
    elif integration_method_lower == 'rk23':
        integration_method_int = 0
    elif integration_method_lower == 'dop853':
        integration_method_int = 2
    else:
        raise UnknownModelError(f"Unsupported integration method provided: {integration_method_lower}.")
    
    cdef str eos_integration_method_lower = eos_integration_method.lower()
    cdef int eos_integration_method_int
    if eos_integration_method_lower == 'rk45':
        eos_integration_method_int = 1
    elif eos_integration_method_lower == 'rk23':
        eos_integration_method_int = 0
    elif eos_integration_method_lower == 'dop853':
        eos_integration_method_int = 2
    else:
        raise UnknownModelError(f"Unsupported EOS integration method provided: {eos_integration_method_lower}.")

    # Convert EOS methods from string to int
    cdef str eos_method_str
    cdef vector[int] eos_integration_method_int_bylayer = vector[int]()
    eos_integration_method_int_bylayer.resize(num_layers)
    cdef int* eos_integration_method_int_bylayer_ptr = eos_integration_method_int_bylayer.data()

    if eos_method_bylayer is None:
        # Use default method
        for layer_i in range(num_layers):
            eos_integration_method_int_bylayer[layer_i] = EOS_INTERPOLATE_METHOD_INT
    else:
        # Step through each layer and check the method
        for layer_i in range(num_layers):
            eos_method_str = eos_method_bylayer[layer_i].lower()
            if eos_method_str == 'interpolate':
                eos_integration_method_int_bylayer[layer_i] = EOS_INTERPOLATE_METHOD_INT
            else:
                raise NotImplementedError("Unknown EOS method provided.")
    
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
            raise ArgumentException(
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
                raise ArgumentException(f"Unsupported value provided for `solve_for`: {solve_for_tmp}.")
    
    # TODO: For now RadialSolver does not support a robust Equation of State method. 
    # The user must provide density, shear, and bulk arrays which will be used to find pressure and gravity.
    cdef double* radius_array_ptr  = &radius_array[0]
    cdef double* density_array_ptr = &density_array[0]

    # Convert complex-valued arrays to C++ complex pointers
    cdef double complex* complex_shear_modulus_ptr = <double complex*> &complex_shear_modulus_array[0]
    cdef double complex* complex_bulk_modulus_ptr  = <double complex*> &complex_bulk_modulus_array[0]

    # Build solution storage
    cdef RadialSolverSolution solution = RadialSolverSolution(
        num_bc_models,
        upper_radius_bylayer_array,
        radius_array,
        degree_l
        )

    # Set the number and type of surface boundary conditions ("tidal", "loading", etc.) that will be solved for 
    # simultaneously.
    solution.set_model_names(bc_models_ptr)

    # Run TidalPy's radial solver function
    cdef rs_error_code = 0
    rs_error_code = cf_radial_solver(
        solution.solution_storage_sptr,
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
    
    # Finalize radial solver solution storage
    solution.finalize_python_storage()

    if log_info:
        solution.print_diagnostics(print_diagnostics = False, log_diagnostics = True)

    if ((not solution.success) or (rs_error_code < 0)) and raise_on_fail:
        if "not implemented" in solution.message:
            raise NotImplementedError(solution.message)
        else:
            raise SolutionFailedError(solution.message)

    if warnings:
        if np.any(solution.steps_taken > 7_000):
            log.warning(f"Large number of steps taken found in radial solver solution (max = {np.max(solution.steps_taken)}). Recommend checking for instabilities (a good method is looking at `<solution>.plot_ys()`).")

    return solution
