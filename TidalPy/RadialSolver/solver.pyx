# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False


from libc.stdio cimport printf
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.math cimport fabs, NAN
from libc.string cimport strcpy

from CyRK cimport PreEvalFunc
from CyRK.utils.vector cimport vector
from CyRK.utils.memory cimport shared_ptr

from TidalPy.logger import get_logger
from TidalPy.exceptions import UnknownModelError

from TidalPy.constants cimport d_G, d_MIN_FREQUENCY, d_MAX_FREQUENCY
from TidalPy.utilities.dimensions.nondimensional cimport (
    cf_non_dimensionalize_physicals,
    cf_redimensionalize_physicals,
    cf_redimensionalize_radial_functions
    )
from TidalPy.RadialSolver.rs_solution cimport RadialSolverSolution
from TidalPy.RadialSolver.shooting cimport cf_shooting_solver
from TidalPy.RadialSolver.matrix cimport cf_matrix_propagate

# EOS Imports (note these will change in a future release)
from TidalPy.Material.eos cimport EOS_ODEInput, solve_eos
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.Material.eos.methods.interpolate cimport InterpolateEOSInput, preeval_interpolate


ctypedef EOS_ODEInput* EOS_ODEInputPtr
ctypedef void* VoidPtr

log = get_logger(__name__)


cdef void cf_radial_solver(
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
        int* is_static_by_layer_ptr,
        int* is_incompressible_by_layer_ptr,
        double surface_pressure,
        unsigned int degree_l,
        size_t num_bc_models,
        int* bc_models_ptr,
        unsigned char core_condition,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        unsigned char integration_method_int,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_by_layer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool nondimensionalize,
        cpp_bool use_prop_matrix,
        unsigned int eos_integration_method,
        double eos_rtol,
        double eos_atol,
        double eos_pressure_tol,
        unsigned int eos_max_iters,
        cpp_bool verbose,
        cpp_bool warnings,
        cpp_bool raise_on_fail,
        ) noexcept nogil:
    
    cdef size_t layer_i, slice_i
    printf("DEBUG-cf_radial_solver Start\n")

    # Feedback
    cdef char[256] message
    cdef char* message_ptr = &message[0]
    
    # Figure out how many slices are in each layer
    cdef vector[size_t] first_slice_index_by_layer_vec = vector[size_t]()
    first_slice_index_by_layer_vec.reserve(num_layers)
    cdef size_t* first_slice_index_by_layer_ptr = &first_slice_index_by_layer_vec[0]

    cdef vector[size_t] num_slices_by_layer_vec = vector[size_t]()
    num_slices_by_layer_vec.reserve(num_layers)
    cdef size_t* num_slices_by_layer_ptr = &num_slices_by_layer_vec[0]

    cdef size_t layer_slices       = 0
    cdef double radius_check       = NAN
    cdef double layer_upper_radius = NAN

    # Pull out raw pointers to avoid repeated calls to the getter
    cdef RadialSolutionStorageCC* solution_storage_ptr = solution_storage_sptr.get()
    cdef EOSSolutionCC* eos_solution_storage_ptr       = solution_storage_sptr.get().eos_solution_sptr.get()

    # Ensure there is at least one layer.
    printf("DEBUG-cf_radial_solver Check on Layers\n")
    if num_layers <= 0:
        solution_storage_ptr.error_code = -5
        strcpy(message_ptr, 'RadialSolver:: requires at least one layer, zero provided.\n')
        solution_storage_ptr.set_message(message_ptr)
        if raise_on_fail:
            printf(message_ptr)
            exit(EXIT_FAILURE)

    printf("DEBUG-cf_radial_solver Start parsing slice index \n")
    if solution_storage_ptr.error_code == 0:
        for layer_i in range(num_layers):
            if layer_i == 0:
                first_slice_index_by_layer_ptr[layer_i] = 0
            else:
                first_slice_index_by_layer_ptr[layer_i] = first_slice_index_by_layer_ptr[layer_i - 1] + 1
            layer_upper_radius = eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i]

            layer_slices = 0
            for slice_i in range(first_slice_index_by_layer_ptr[layer_i], total_slices):
                radius_check = radius_array_in_ptr[slice_i]
                if radius_check > layer_upper_radius:
                    # We have passed this layer.
                    break
                else:
                    layer_slices += 1
            
            if layer_slices <= 3:
                solution_storage_ptr.error_code == -5
                strcpy(message_ptr, 'RadialSolver:: At least three layer slices per layer are required. Try using more slices in the input arrays.\n')
                solution_storage_sptr.get().set_message(message_ptr)
                if verbose or raise_on_fail:
                    printf(message_ptr)
                if raise_on_fail:
                    exit(EXIT_FAILURE)

            num_slices_by_layer_ptr[layer_i] = layer_slices

    # Get other needed inputs
    cdef double radius_planet = radius_array_in_ptr[total_slices - 1]

    # Non-dimensionalize inputs
    cdef double G_to_use                = NAN
    cdef double radius_planet_to_use    = NAN
    cdef double bulk_density_to_use     = NAN
    cdef double frequency_to_use        = NAN
    cdef double surface_pressure_to_use = NAN
    
    if nondimensionalize and solution_storage_ptr.error_code == 0:
        cf_non_dimensionalize_physicals(
            total_slices,
            frequency,
            radius_planet,
            planet_bulk_density,
            radius_array_in_ptr,
            density_array_in_ptr,
            NULL,  # Pressure ptr (not found yet.)
            NULL,  # Gravity ptr (not found yet.)
            complex_bulk_modulus_in_ptr,
            complex_shear_modulus_in_ptr,
            &radius_planet_to_use,
            &bulk_density_to_use,
            &surface_pressure_to_use,
            &frequency_to_use,
            &G_to_use
            )

        for layer_i in range(num_layers):
            eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i] /= radius_planet
        
        # Change starting radius if it was provided (if it is set to default then it is 0 and this won't have an effect)
        starting_radius /= radius_planet

        # Update the radius array inside the C++ classes
        solution_storage_ptr.change_radius_array(radius_array_in_ptr, total_slices, True)

    else:
        G_to_use                = d_G
        radius_planet_to_use    = radius_planet
        bulk_density_to_use     = planet_bulk_density
        frequency_to_use        = frequency
        surface_pressure_to_use = surface_pressure

    # Solve the equaiton of state for the planet

    # TODO: For now there is only one accepted EOS, the interpolated kind. In the future additional EOS will be supplied
    # either via arguments to this function or a more OOP approach where they are built into the layers.
    # Build arrays of EOS inputs.
    printf("DEBUG-cf_radial_solver EOS setup\n")
    cdef size_t bottom_slice_index
    cdef vector[PreEvalFunc] eos_function_bylayer_vec = vector[PreEvalFunc]()
    eos_function_bylayer_vec.reserve(num_layers)

    # Build vector of EOS inputs
    cdef EOS_ODEInput eos_input
    cdef vector[EOS_ODEInput] eos_inputs_bylayer_vec = vector[EOS_ODEInput]()
    eos_inputs_bylayer_vec.reserve(num_layers)

    # TODO: Below is specific to interpolate EOS
    cdef vector[void] specific_eos_input_bylayer_vec = vector[void]()
    specific_eos_input_bylayer_vec.reserve(num_layers)

    # Build Equation of State functions and input data structures. Record memory addresses for use by the EOS solver
    if solution_storage_ptr.error_code == 0:
        for layer_i in range(num_layers):
            # TODO: Below is specific to interpolate EOS. For now we are only storing the interpolate version of the EOS for each layer.
            eos_function_bylayer_vec.push_back(preeval_interpolate)

            # Build EOS input
            bottom_slice_index = first_slice_index_by_layer_ptr[layer_i]

            # Build EOS input for specific EOS model
            # TODO: Below is specific to interpolate EOS.
            specific_eos_input_bylayer_vec.emplace_back(
                num_slices_by_layer_vec[layer_i],                   # Number of slices for this layer [size_t]
                &radius_array_in_ptr[bottom_slice_index],           # Radius array pointer [double*]
                &density_array_in_ptr[bottom_slice_index],          # Density array pointer [double*]
                &complex_bulk_modulus_in_ptr[bottom_slice_index],   # Complex bulk array pointer [double complex*]
                &complex_shear_modulus_in_ptr[bottom_slice_index],  # Complex shear array pointer [double complex*]
                )
            
            # Build input for generalized EOS solver
            eos_inputs_bylayer_vec.emplace_back(
                G_to_use,                                        # Gravitational constant [double]
                radius_planet_to_use,                            # Planet radius [double]
                specific_eos_input_bylayer_vec[layer_i],         # void-casted pointer to model specific input (created just above)
                False,                                           # Final solve flag [bool] (will be updated by EOS solver)
                False,                                           # Final update shear flag [bool] (will be updated by EOS solver)
                False                                            # Final update bulk flag [bool] (will be updated by EOS solver)
            )
            # TODO: update bulk/shear flags are overwritten by EOS solver regardless of layer type. They should not ever need to be updated, for example shear for liquid layers?
        
    if solution_storage_ptr.error_code == 0:
        printf("DEBUG-cf_radial_solver EOS Solver\n")
        solve_eos(
            solution_storage_ptr.eos_solution_sptr,  # Equation of state storage C++ class
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
    
    printf("DEBUG::radial_solver - EOS Results:\n")
    printf("\t EOS Solution Sptr            = %p\n", solution_storage_ptr.eos_solution_sptr)
    printf("\t EOS Solution ptr             = %p\n", solution_storage_ptr.eos_solution_sptr.get())
    printf("\t EOS CySolver sptr            = %p\n", solution_storage_ptr.eos_solution_sptr.get().cysolver_results_sptr_bylayer_vec[0])
    printf("\t EOS CySolver ptr (from sptr) = %p\n", solution_storage_ptr.eos_solution_sptr.get().cysolver_results_sptr_bylayer_vec[0].get())
    printf("\t EOS CySolver ptr (stored)    = %p\n", solution_storage_ptr.eos_solution_sptr.get().cysolver_results_ptr_bylayer_vec[0])
    printf("\t\t mass = %e\n", solution_storage_ptr.eos_solution_sptr.get().mass)
    printf("\t\t moi  = %e\n", solution_storage_ptr.eos_solution_sptr.get().moi)
    printf("\t\t g    = %e\n", solution_storage_ptr.eos_solution_sptr.get().surface_gravity)
    printf("\t\t Psur = %e\n", solution_storage_ptr.eos_solution_sptr.get().surface_pressure)
    printf("\t\t R    = %e\n", solution_storage_ptr.eos_solution_sptr.get().radius)
    printf("\t\t Pcen = %e\n", solution_storage_ptr.eos_solution_sptr.get().central_pressure)

    exit(-86)

    # Step through the radial steps to find EOS-dependent parameters
    if eos_solution_storage_ptr.success and solution_storage_ptr.error_code == 0:
        # Run requested radial solver method
        if use_prop_matrix:

            cf_matrix_propagate(
                solution_storage_sptr,     # (Modified) Final radial solution storage struct pointer [RadialSolutionStorageCC*]
                frequency_to_use,          # Forcing frequency [double]
                bulk_density_to_use,       # Planet bulk density [double]
                # TODO: In the future the propagation matrix should take in layer types and multiple layers
                # int* layer_types_ptr,
                # int* is_static_by_layer_ptr,
                # int* is_incompressible_by_layer_ptr,
                num_bc_models,             # Number of boundary conditions requested by user [size_t]
                bc_models_ptr,             # Boundary condition model int array pointer [int*]
                G_to_use,                  # Gravitational constant [double]
                degree_l,                  # Harmonic degree [unsigned int]
                core_condition,            # Starting condition model int at the inner boundary (usually a core) see TidalPy.RadialSolver.matrix.pyx for options [unsigned char]
                verbose,                   # Verbose flag [cpp_bool]
                raise_on_fail              # Flag to allow for early crashes when integration fails [cpp_bool]
                )
        else:
            printf("DEBUG-cf_radial_solver - Pre shooting method\n")
            cf_shooting_solver(
                solution_storage_sptr,          # (Modified) Final radial solution storage struct pointer [RadialSolutionStorageCC*]
                frequency_to_use,               # Forcing frequency [double]
                bulk_density_to_use,            # Planet bulk density [double]
                layer_types_ptr,                # Layer type int  array pointer [int*]
                is_static_by_layer_ptr,         # Layer is_static flag array pointer [int*]
                is_incompressible_by_layer_ptr, # Pointer array of layer is_incompressible flag array pointer [int*]
                first_slice_index_by_layer_ptr, # First radial slice of each layer array pointer [size_t*]
                num_slices_by_layer_ptr,        # Number of radial slices in each layer array pointer [size_t*]
                num_bc_models,                  # Number of boundary conditions requested by user [size_t]
                bc_models_ptr,                  # Boundary condition model int array pointer [int*]
                G_to_use,                       # Gravitational constant [double]
                degree_l,                       # Harmonic degree [unsigned int]
                use_kamata,                     # Flag to use Kamata+ (2015)'s starting conditions vs. Takeuchi+Saito (1972) [cpp_bool]
                starting_radius,                # Starting radius for solver. For higher degree solutions you generally want to start higher up in the planet. [double]
                start_radius_tolerance,         # Tolerance used if `starting_radius` is not provided. [double]
                integration_method_int,         # Integration method int (0=RK23, 1=RK45, 2=DOP853) [unsigned char]
                integration_rtol,               # Integration relative tolerance [double]
                integration_atol,               # Integration absolute tolerance [double]
                scale_rtols_by_layer_type,      # Flag for if tolerances should vary with layer type (using pre-defined scaling) [cpp_bool]
                max_num_steps,                  # Maximum number of integration steps allowed [size_t]
                expected_size,                  # Expected number of integration steps required for the average layer [size_t]
                max_ram_MB,                     # Maximum amount of ram allowed for each layer's integration (note if parallized then radial solver will exceed this value; there is also overhead of other functions) [size_t]
                max_step,                       # Maximum allowed step size per layer [double]
                verbose,                        # Verbose flag [cpp_bool]
                raise_on_fail                   # Flag to allow for early crashes when integration fails [cpp_bool]
                )
            printf("DEBUG-cf_radial_solver - Post shooting method\n")

        # FIX ME: left off
    # Finalize solution storage
    # Get a reference pointer to solution array
    cdef double* solution_dbl_ptr = solution_storage_ptr.full_solution_ptr
    # Cast the solution pointer from double to double complex
    cdef double complex* solution_ptr = <double complex*>solution_dbl_ptr

    if nondimensionalize:
        # Redimensionalize eos properties
        cf_redimensionalize_physicals(
            total_slices,
            frequency,
            radius_planet,
            planet_bulk_density,
            solution_storage_ptr.radius_array_ptr,
            solution_storage_ptr.density_array_ptr,
            solution_storage_ptr.pressure_array_ptr,
            solution_storage_ptr.gravity_array_ptr,
            solution_storage_ptr.complex_bulk_array_ptr,
            solution_storage_ptr.complex_shear_array_ptr,
            &radius_planet_to_use,
            &bulk_density_to_use,
            &frequency_to_use,
            &G_to_use
            )

        for layer_i in range(num_layers):
            eos_solution_storage_ptr.upper_radius_bylayer_vec[layer_i] *= radius_planet

    if solution_storage_ptr.success:
        printf("DEBUG-cf_radial_solver - Successful closeout\n")
        if nondimensionalize:
            # Redimensionalize the solution 
            cf_redimensionalize_radial_functions(
                solution_ptr,
                radius_planet,
                planet_bulk_density,
                total_slices,
                num_bc_models)

        # Calculate Love numbers
        printf("DEBUG-cf_radial_solver - Pre find-love\n")
        solution_storage_ptr.find_love()
        printf("DEBUG-cf_radial_solver - Post find-love\n")


def radial_solver(
        double[::1] radius_array,
        double[::1] density_array,
        double complex[::1] complex_bulk_modulus_array,
        double complex[::1] complex_shear_modulus_array,
        double frequency,
        double planet_bulk_density,
        tuple layer_types,
        tuple is_static_by_layer,
        tuple is_incompressible_by_layer,
        double[::1] upper_radius_by_layer_array,
        double surface_pressure = 0.0,
        unsigned int degree_l = 2,
        tuple solve_for = None,
        unsigned char core_condition = 0,
        cpp_bool use_kamata = False,
        double starting_radius = 0.0,
        double start_radius_tolerance = 1.0e-5,
        str integration_method = 'RK45',
        double integration_rtol = 1.0e-6,
        double integration_atol = 1.0e-12,
        cpp_bool scale_rtols_by_layer_type = False,
        size_t max_num_steps = 500_000,
        size_t expected_size = 500,
        size_t max_ram_MB = 500,
        double max_step = 0,
        cpp_bool nondimensionalize = True,
        cpp_bool use_prop_matrix = False,
        str eos_integration_method = 'RK45',
        double eos_rtol = 1.0e-3,
        double eos_atol = 1.0e-5,
        double eos_pressure_tol = 1.0e-3,
        unsigned int eos_max_iters = 40,
        cpp_bool verbose = False,
        cpp_bool warnings = True,
        cpp_bool raise_on_fail = False,
        cpp_bool perform_checks = True
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
    verbose : bool, default=False
        If True, then additioal information will be printed to the terminal during the solution. 
    warnings : bool, default=True
        If True, then warnings will be printed to the terminal during the solution. 
    raise_on_fail : bool, default=False
        If Ture, then the solver will raise an exception if integration was not successful. By default RadialSolver
        fails silently.
    perform_checks : bool, default=True
        Performs sanity checks that raise python exceptions. If turned off then these checks will be skipped providing 
        some boost to performance but at the risk of uncaught exceptions (crashes).
    
    Returns
    -------
    solution : RadialSolverSolution
        Solution to the viscoelastic-gravitational problem inside a planet.
        Also contains the EOS solver solution for the entire planet.
    """

    printf("DEBUG-RadialSolver Point 1\n")

    cdef size_t total_slices = radius_array.size

    # Unpack inefficient user-provided tuples into bool arrays and pass by pointer
    cdef size_t num_layers = len(layer_types)

    # Check on solver type
    if use_prop_matrix and num_layers > 1:
        raise NotImplementedError("Currently, TidalPy's propagation matrix technique only works for 1-layer worlds. For 2 layer worlds where the lower layer is a liquid: you can start the solver at the bottom of the upper solid layer.")

    printf("DEBUG-RadialSolver Point 2\n")
    if perform_checks:
        assert density_array.size               == total_slices
        assert complex_bulk_modulus_array.size  == total_slices
        assert complex_shear_modulus_array.size == total_slices
        # Check that number of assumptions match.
        if len(is_static_by_layer) != num_layers:
            raise AttributeError('Number of `is_static_by_layer` must match number of `layer_types`.')
        if len(is_incompressible_by_layer) != num_layers:
            raise AttributeError('Number of `is_incompressible_by_layer` must match number of `layer_types`.')
        if upper_radius_by_layer_array.size != num_layers:
            raise AttributeError('Number of `upper_radius_by_layer` must match number of `layer_types`.')
        if radius_array[0] != 0.:
            raise AttributeError('Radius array must start at zero.')

        if fabs(frequency) < d_MIN_FREQUENCY:
            raise ValueError('Forcing frequency is too small (are you sure you are in rad s-1?).')
        elif fabs(frequency) > d_MAX_FREQUENCY:
            raise ValueError('Forcing frequency is too large (are you sure you are in rad s-1?).')
        
    printf("DEBUG-RadialSolver Point 3\n")
    # Build array of assumptions
    # OPT: Perhaps set a maximum number of layers then we can put these on the stack rather than heap allocating them.
    cdef vector[int] layer_assumptions_vec   = vector[int]()
    layer_assumptions_vec.reserve(3 * num_layers)
    cdef int* layer_assumptions_ptr          = &layer_assumptions_vec[0]

    cdef int* layer_types_ptr                = &layer_assumptions_ptr[0]
    cdef int* is_static_by_layer_ptr         = &layer_assumptions_ptr[num_layers]
    cdef int* is_incompressible_by_layer_ptr = &layer_assumptions_ptr[2 * num_layers]
    printf("DEBUG-RadialSolver Point 4\n")

    cdef str layer_type
    cdef cpp_bool dynamic_liquid = False

    # Pull out information for each layer and store in heap memory
    printf("DEBUG-RadialSolver Point 6\n")
    cdef size_t layer_i
    for layer_i in range(num_layers):
        layer_type = layer_types[layer_i]
        
        is_static_by_layer_ptr[layer_i]         = is_static_by_layer[layer_i]
        is_incompressible_by_layer_ptr[layer_i] = is_incompressible_by_layer[layer_i]

        if not dynamic_liquid:
            if (layer_type == 1) and not is_static_by_layer_ptr[layer_i]:
                # There is at least one dynamic liquid layer
                dynamic_liquid = True

        # Convert user-provided strings to ints for the layer type
        if layer_type.lower() == 'solid':
            layer_types_ptr[layer_i] = 0
        elif layer_type.lower() == 'liquid':
            layer_types_ptr[layer_i] = 1
        else:
            layer_types_ptr[layer_i] = -1
            log.error(f"Layer type {layer_type} is not supported. Currently supported types: 'solid', 'liquid'.")
            raise UnknownModelError(f"Layer type {layer_type} is not supported. Currently supported types: 'solid', 'liquid'.")
    
    printf("DEBUG-RadialSolver Point 7\n")
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
    printf("DEBUG-RadialSolver Point 8\n")
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
    
    cdef str eos_integration_method_lower = eos_integration_method.lower()
    cdef unsigned char eos_integration_method_int
    if eos_integration_method_lower == 'rk45':
        eos_integration_method_int = 1
    elif eos_integration_method_lower == 'rk23':
        eos_integration_method_int = 0
    elif eos_integration_method_lower == 'dop853':
        eos_integration_method_int = 2
    else:
        log.error(f"Unsupported EOS integration method provided: {eos_integration_method_lower}.")
        raise UnknownModelError(f"Unsupported EOS integration method provided: {eos_integration_method_lower}.")
    
    # Clean up what values the solver is solving for.
    printf("DEBUG-RadialSolver Point 9\n")
    cdef int[5] bc_models
    cdef size_t num_bc_models
    cdef int* bc_models_ptr = &bc_models[0]
    cdef str solve_for_tmp 

    printf("DEBUG-RadialSolver Point 10\n")
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
    
    # TODO: For now RadialSolver does not support a robust Equation of State method. 
    # The user must provide density, shear, and bulk arrays which will be used to find pressure and gravity.
    cdef double* radius_array_ptr  = &radius_array[0]
    cdef double* density_array_ptr = &density_array[0]

    # Convert complex-valued arrays to C++ complex pointers
    cdef double complex* complex_shear_modulus_ptr = <double complex*> &complex_shear_modulus_array[0]
    cdef double complex* complex_bulk_modulus_ptr  = <double complex*> &complex_bulk_modulus_array[0]

    printf("DEBUG-RadialSolver Point 11\n")
    # Build solution storage
    cdef RadialSolverSolution solution = RadialSolverSolution(num_bc_models, upper_radius_by_layer_array, radius_array)
    
    printf("DEBUG-RadialSolver Point 11b\n")
    solution.set_model_names(bc_models_ptr)
    printf("DEBUG-RadialSolver Point 11c\n")

    # Run TidalPy's radial solver function
    printf("DEBUG-RadialSolver Point 13a - Pre cf_radial_solver call\n")
    cf_radial_solver(
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
        is_static_by_layer_ptr,
        is_incompressible_by_layer_ptr,
        surface_pressure,
        degree_l,
        num_bc_models,
        bc_models_ptr,
        core_condition,
        use_kamata,
        starting_radius,
        start_radius_tolerance,
        integration_method_int,
        integration_rtol,
        integration_atol,
        scale_rtols_by_layer_type,
        max_num_steps,
        expected_size,
        max_ram_MB,
        max_step,
        nondimensionalize,
        use_prop_matrix,
        eos_integration_method_int,
        eos_rtol,
        eos_atol,
        eos_pressure_tol,
        eos_max_iters,
        verbose,
        warnings,
        raise_on_fail,
        )
    printf("DEBUG-RadialSolver Point 13b - Post cf_radial_solver call\n")

    return solution
