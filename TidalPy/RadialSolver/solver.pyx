# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as np

from TidalPy.logger import get_logger
from TidalPy.exceptions import UnknownModelError
log = get_logger(__name__)

from libc.math cimport fabs, NAN
from CyRK.utils.utils cimport allocate_mem, reallocate_mem, free_mem
from TidalPy.RadialSolver.shooting cimport cf_shooting_solver

# Maximum size for array building
cdef size_t MAX_NUM_Y      = 6
cdef size_t MAX_NUM_Y_REAL = 2 * MAX_NUM_Y
cdef size_t MAX_NUM_SOL    = 3


cdef class RadialSolverSolution():

    def __init__(
            self,
            size_t num_slices,
            size_t* bc_models_ptr,
            size_t num_ytypes
            ):

        # loop indicies
        cdef size_t i
        cdef size_t love_array_size

        # Initialize pointers
        self.full_solution_ptr = NULL
        self.complex_love_ptr = NULL

        # Initialize status
        self._message = 'RadialSolverSolution has not had its status set.'
        self.success = False

        # Store number of solution types
        cdef size_t bc_model
        self.num_ytypes = num_ytypes

        for i in range(self.num_ytypes):
            bc_model = bc_models_ptr[i]
            if bc_model == 0:
                self.ytypes[i] = "free"
            elif bc_model == 1:
                self.ytypes[i] = "tidal"
            elif bc_model == 2:
                self.ytypes[i] = "loading"
            else:
                raise AttributeError(f"Unknown boundary condition {bc_model}")
        
        # Store size information
        self.num_slices = num_slices
        self.total_size = MAX_NUM_Y * self.num_slices * self.num_ytypes

        # Have the radial solver take control of the full solution memory
        self.full_solution_ptr = <double complex*> allocate_mem(
            self.total_size * sizeof(double complex),
            'full_solution_ptr (RadialSolverSolution; init)'
            )
        if not (self.full_solution_ptr is NULL):
            self.full_solution_view = <double complex[:self.total_size]> self.full_solution_ptr
        
        # Set all values of the solution array to NANs. This ensures that if there is a problem with the solver then
        #  the solutions will be defined (but nan).
        for i in range(self.total_size):
            self.full_solution_ptr[i] = NAN
        
        # Initialize Love numbers (3 love numbers for each requested y-type)
        # Love numbers are stored (k, h, l)_ytype0, (k, h, l)_ytype1, (k, h, l)_ytype2, ...
        love_array_size = 3 * self.num_ytypes
        self.complex_love_ptr = <double complex*> allocate_mem(
            love_array_size * sizeof(double complex),
            'complex_love_ptr (RadialSolverSolution; init)'
            )
        if not (self.complex_love_ptr is NULL):
            self.complex_love_view = <double complex[:love_array_size]> self.complex_love_ptr
        for i in range(love_array_size):
            self.complex_love_ptr[i] = NAN

    @property
    def message(self):
        """ Return solver's message """
        return str(self._message, 'UTF-8')
    
    @property
    def result(self):
        """ Return result array. """

        if self.success:
            return np.ascontiguousarray(
                self.full_solution_view,
                dtype=np.complex128
                ).reshape((self.num_slices, self.num_ytypes * MAX_NUM_Y)).T
        else:
            return None
        
    @property
    def love(self):
        """ Return all complex love numbers. """
        if self.success:
            return np.ascontiguousarray(
                self.complex_love_view,
                dtype=np.complex128
            ).reshape((self.num_ytypes, 3))
        else:
            return None

    @property
    def k(self):
        """ Tidal Love number k. """
        if self.success:
            return np.ascontiguousarray(
                self.complex_love_view[0::3],
                dtype=np.complex128
            )
        else:
            return None

    @property
    def h(self):
        """ Tidal Love number h. """
        if self.success:
            return np.ascontiguousarray(
                self.complex_love_view[1::3],
                dtype=np.complex128
            )
        else:
            return None
    
    @property
    def l(self):
        """ Tidal Shida number l. """
        if self.success:
            return np.ascontiguousarray(
                self.complex_love_view[2::3],
                dtype=np.complex128
            )
        else:
            return None
    
    def __len__(self):
        """Return number of solution types."""
        return <Py_ssize_t>self.num_ytypes
    
    def __getitem__(self, str ytype_name):
        """Get a specific solution type array."""
        
        cdef size_t i
        cdef size_t requested_sol_num = 0
        cdef bint found = False
        cdef str sol_test_name
        for i in range(self.num_ytypes):
            sol_test_name = str(self.ytypes[i], 'UTF-8')
            if sol_test_name == ytype_name:
                requested_sol_num = i
                found = True
                break
        if not found:
            raise AttributeError('Unknown solution type requested. Key must match names passed to radial_solver "solve_for" argument and be lower case.')
        
        # Slice the result and return only the requested solution type.
        if self.success:
            return self.result[MAX_NUM_Y * (requested_sol_num): MAX_NUM_Y * (requested_sol_num + 1)]
        else:
            return None

    def __dealloc__(self):

        # The RadialSolverSolution class has full control of the solution so it is responsible for releasing its memory.
        if not (self.full_solution_ptr is NULL):
            free_mem(self.full_solution_ptr)
        if not (self.complex_love_ptr is NULL):
            free_mem(self.complex_love_ptr)


def radial_solver(
        double[::1] radius_array,
        double[::1] density_array,
        double[::1] gravity_array,
        double[::1] bulk_modulus_array,
        double complex[::1] complex_shear_modulus_array,
        double frequency,
        double planet_bulk_density,
        tuple layer_types,
        tuple is_static_by_layer,
        tuple is_incompressible_by_layer,
        tuple upper_radius_by_layer,
        unsigned int degree_l = 2,
        tuple solve_for = None,
        bint use_kamata = False,
        str integration_method = 'RK45',
        double integration_rtol = 1.0e-6,
        double integration_atol = 1.0e-12,
        bint scale_rtols_by_layer_type = False,
        size_t max_num_steps = 500_000,
        size_t expected_size = 500,
        size_t max_ram_MB = 500,
        double max_step = 0,
        bint limit_solution_to_radius = True,
        bint nondimensionalize = True,
        bint use_prop_matrix = False,
        bint verbose = False,
        bint warnings = True,
        bint raise_on_fail = False
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
    gravity_array : np.ndarray[dtype=np.float64]
        Acceleration due to gravity at each radius [m-2].
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
    assert gravity_array.size               == total_slices
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
    cdef size_t[5] bc_models
    cdef size_t num_bc_models
    cdef size_t* bc_models_ptr = &bc_models[0]
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

    # Prepare to run
    cdef RadialSolverSolution result
    try:
        result = cf_shooting_solver(
            total_slices,
            &radius_array[0],
            &density_array[0],
            &gravity_array[0],
            &bulk_modulus_array[0],
            &complex_shear_modulus_array[0],
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
            raise_on_fail)
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
    
    return result
