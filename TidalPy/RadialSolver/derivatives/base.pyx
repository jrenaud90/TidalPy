# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from libc.math cimport pi

from CyRK.cy.common cimport MAX_STEP, EPS_100
from CyRK.cy.cysolver cimport CySolver
from CyRK.array.interp cimport interpj_ptr, interp_ptr, interp_complex_ptr

cdef class RadialSolverBase(CySolver):
    def __init__(
            self,
            # RadialSolverBase Inputs
            double frequency,
            unsigned int degree_l,
            double G_to_use,

            # Regular CySolver Inputs
            const double[::1] y0,
            (double, double) t_span,
            unsigned char rk_method = 1,
            double max_step = MAX_STEP,
            double first_step = 0.,
            size_t max_num_steps = 0,
            size_t expected_size = 0,
            size_t max_ram_MB = 500
            ):

        # Initialize pointers to null
        self.radius_array_ptr = NULL
        self.density_array_ptr = NULL
        self.gravity_array_ptr = NULL
        self.bulk_modulus_array_ptr = NULL
        self.shear_modulus_array_ptr = NULL

        # Load in floats and ints
        self.frequency  = frequency
        self.degree_l   = degree_l
        self.G_to_use   = G_to_use
        self.grav_coeff = 4. * pi * self.G_to_use

        # Initialize state variables
        self.shear_modulus = 0. + 0.j
        self.bulk_modulus  = 0.
        self.density = 0.
        self.gravity = 0.

        # Setup regular CySolver
        super().__init__(
            t_span=t_span,
            y0=y0,
            args=None,
            rk_method=rk_method,
            max_step=max_step,
            first_step=first_step,
            max_num_steps=max_num_steps,
            t_eval=None,
            capture_extra=False,
            num_extra=0,
            interpolate_extra=False,
            expected_size=expected_size,
            max_ram_MB=max_ram_MB,
            call_first_reset=False,
            auto_solve=False
            )


    cdef void install_pointers(
            self,

            # RadialSolverBase pointers
            size_t num_slices,
            double* radius_array_ptr,
            double* density_array_ptr,
            double* gravity_array_ptr,
            double* bulk_modulus_array_ptr,
            double complex* shear_modulus_array_ptr,
            double* atols,
            double* rtols,

            # Additional optional arguments for RadialSolver class
            bint limit_solution_to_radius = True,
            bint call_first_reset = False,
            bint auto_solve = True,
            ):
        # Cython does not support non-python objects being passed to __init__ or __cinit__. So we need this helper
        # method to take the required pointers and load them into the class.

        # Setup loop variables
        cdef size_t i

        self.num_slices = num_slices
        # Set class pointers to inputs
        self.radius_array_ptr        = radius_array_ptr
        self.density_array_ptr       = density_array_ptr
        self.gravity_array_ptr       = gravity_array_ptr
        self.bulk_modulus_array_ptr  = bulk_modulus_array_ptr
        self.shear_modulus_array_ptr = shear_modulus_array_ptr

        # Determine if solution should only be provided at points in the provided radius array.
        if limit_solution_to_radius:
            # The t_eval pointer is allocated by the CySolver class (even if t_eval=None) was passed during init.
            # We want to use t_eval = radius_ptr but we can not just set it like that because:
            #   1) we will lose the reference to the original t_eval_ptr that was allocated, leading to a memory leak.
            #   2) when the solver class is dealloc it will dealloc the radius_ptr which we don't want the solver to own.
            # Instead we will realloc the memory and set its values equal to radius pointer.
            self.change_t_eval_pointer(self.radius_array_ptr, num_slices, auto_reset_state=False)

        # rtols and atols are provided as pointers which the base class does not support right away.
        # Set those up now.
        cdef double rtol_tmp
        for i in range(self.y_size):
            rtol_tmp = rtols[i]
            if rtol_tmp < EPS_100:
                rtol_tmp = EPS_100
            self.rtols_ptr[i] = rtol_tmp
            self.atols_ptr[i] = atols[i]

        # Reset the state (this will be the first time since we passed "call_first_reset=False" to the parent class)
        if call_first_reset or auto_solve:
            self.reset_state()

        # Run integrator if requested.
        if auto_solve:
            self._solve(reset=False)

    cdef void update_constants(self) noexcept nogil:

        # Define degree l constants.
        cdef double degree_l_flt

        degree_l_flt = <double>self.degree_l
        self.lp1     = degree_l_flt + 1.
        self.lm1     = degree_l_flt - 1.
        self.llp1    = degree_l_flt * self.lp1


    cdef void update_interp(
            self,
            bint update_bulk,
            bint update_shear
            ) noexcept nogil:

        # Set state variables based on an interpolation using the provided radius.
        cdef (double, Py_ssize_t) interp_out
        cdef Py_ssize_t index_j

        # The first interpolation will be the slowest as it must find the closest index.
        # We will use this index in the other interpolations.
        interp_out = interpj_ptr(
            self.t_now,
            self.radius_array_ptr,
            self.density_array_ptr,
            self.num_slices,
            provided_j=-2
            )

        self.density = interp_out[0]
        index_j      = interp_out[1]

        self.gravity = interp_ptr(
            self.t_now,
            self.radius_array_ptr,
            self.gravity_array_ptr,
            self.num_slices,
            provided_j=index_j
            )

        if update_bulk:
            self.bulk_modulus = interp_ptr(
                self.t_now,
                self.radius_array_ptr,
                self.bulk_modulus_array_ptr,
                self.num_slices,
                provided_j=index_j
                )

        if update_shear:
            self.shear_modulus = interp_complex_ptr(
                self.t_now,
                self.radius_array_ptr,
                self.shear_modulus_array_ptr,
                self.num_slices,
                provided_j=index_j
                )

    def __dealloc__(self):
        # Release memory held by the parent class.
        # Note: we do not need to call the super class' dealloc. Per Cython documentation:
        # "When subclassing extension types, be aware that the __dealloc__() method of the superclass will always be
        #  called, even if it is overridden. This is in contrast to typical Python behavior where superclass methods
        #  will not be executed unless they are explicitly called by the subclass."
        # super().__dealloc__()

        # Release memory held by this class.
        pass