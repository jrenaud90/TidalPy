# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from libcpp cimport bool as bool_cpp_t
from libc.math cimport pi

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from CyRK.cy.cysolver cimport CySolver, MAX_STEP
from CyRK.array.interp cimport interpj_ptr, interp_ptr, interp_complex_ptr


cdef class RadialSolverBase(CySolver):
    def __init__(
            self,

            # RadialSolverBase Inputs
            const double[::1] radius_array_view,
            const double[::1] density_array_view,
            const double[::1] gravity_array_view,
            const double complex[::1] shear_modulus_array_view,
            const double[::1] bulk_modulus_array_view,
            double frequency,
            unsigned int degree_l,
            double G_to_use,

            # Regular CySolver Inputs
            (double, double) t_span,
            const double[::1] y0,
            tuple args = None,
            double rtol = 1.e-3,
            double atol = 1.e-6,
            double[::1] rtols = None,
            double[::1] atols = None,
            unsigned char rk_method = 1,
            double max_step = MAX_STEP,
            double first_step = 0.,
            Py_ssize_t max_num_steps = 0,
            const double[::1] t_eval = None,
            bool_cpp_t capture_extra = False,
            Py_ssize_t num_extra = 0,
            bool_cpp_t interpolate_extra = False,
            Py_ssize_t expected_size = 0,
            bool_cpp_t auto_solve = True,

            # Additional optional arguments for RadialSolver class
            bool_cpp_t limit_solution_to_radius = True
            ):

        # Load in floats and ints
        self.n_radius = len(radius_array_view)

        self.frequency  = frequency
        self.degree_l   = degree_l
        self.G_to_use   = G_to_use
        self.grav_coeff = 4. * pi * self.G_to_use

        # Determine if solution should only be provided at points in the provided radius array.
        if limit_solution_to_radius:
            t_eval = radius_array_view
        else:
            t_eval = None

        # Claim memory for pointer arrays
        cdef Py_ssize_t i
        self.radius_array_ptr = <double *> PyMem_Malloc(self.n_radius * sizeof(double))
        if not self.radius_array_ptr:
            raise MemoryError()

        self.density_array_ptr = <double *> PyMem_Malloc(self.n_radius * sizeof(double))
        if not self.density_array_ptr:
            raise MemoryError()

        self.gravity_array_ptr = <double *> PyMem_Malloc(self.n_radius * sizeof(double))
        if not self.gravity_array_ptr:
            raise MemoryError()

        self.shear_modulus_array_ptr = <double complex*> PyMem_Malloc(self.n_radius * sizeof(double complex))
        if not self.shear_modulus_array_ptr:
            raise MemoryError()

        self.bulk_modulus_array_ptr = <double *> PyMem_Malloc(self.n_radius * sizeof(double))
        if not self.bulk_modulus_array_ptr:
            raise MemoryError()

        # Populate values
        for i in range(self.n_radius):
            self.radius_array_ptr[i]        = radius_array_view[i]
            self.density_array_ptr[i]       = density_array_view[i]
            self.gravity_array_ptr[i]       = gravity_array_view[i]
            self.shear_modulus_array_ptr[i] = shear_modulus_array_view[i]
            self.bulk_modulus_array_ptr[i]  = bulk_modulus_array_view[i]

        # Initialize state variables
        self.shear_modulus = 0. + 0.j
        self.bulk_modulus  = 0.
        self.density = 0.
        self.gravity = 0.

        # Setup regular CySolver
        super().__init__(
            t_span=t_span, y0=y0, args=args,
            rtol=rtol, atol=atol, rtols=rtols, atols=atols, rk_method=rk_method,
            max_step=max_step, first_step=first_step, max_num_steps=max_num_steps,
            t_eval=t_eval,
            capture_extra=capture_extra, num_extra=num_extra, interpolate_extra=interpolate_extra,
            expected_size=expected_size, auto_solve=auto_solve
            )

    cdef void update_constants(self) noexcept nogil:

        # Define degree l constants.
        cdef double degree_l_flt

        degree_l_flt = <double>self.degree_l
        self.lp1     = degree_l_flt + 1.
        self.lm1     = degree_l_flt - 1.
        self.llp1    = degree_l_flt * self.lp1

    cdef void update_interp(
            self,
            bool_cpp_t update_bulk,
            bool_cpp_t update_shear
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
            self.n_radius,
            provided_j=-2
            )

        self.density = interp_out[0]
        index_j      = interp_out[1]

        self.gravity = interp_ptr(
            self.t_now,
            self.radius_array_ptr,
            self.gravity_array_ptr,
            self.n_radius,
            provided_j=index_j
            )

        if update_bulk:
            self.bulk_modulus = interp_ptr(
                self.t_now,
                self.radius_array_ptr,
                self.bulk_modulus_array_ptr,
                self.n_radius,
                provided_j=index_j
                )

        if update_shear:
            self.shear_modulus = interp_complex_ptr(
                self.t_now,
                self.radius_array_ptr,
                self.shear_modulus_array_ptr,
                self.n_radius,
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
        PyMem_Free(self.radius_array_ptr)
        PyMem_Free(self.density_array_ptr)
        PyMem_Free(self.gravity_array_ptr)
        PyMem_Free(self.shear_modulus_array_ptr)
        PyMem_Free(self.bulk_modulus_array_ptr)
