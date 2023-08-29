# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as np
np.import_array()

from libcpp cimport bool as bool_cpp_t
from libc.math cimport pi

from CyRK.cy.cysolver cimport CySolver, MAX_STEP
from CyRK.array.interp cimport interpj, interp, interp_complex

cdef double EPS = np.finfo(np.float64).eps
cdef double EPS_10 = EPS * 10.
cdef double EPS_100 = EPS * 100.


cdef class RadialSolverBase(CySolver):
    def __init__(
            self,

            # RadialSolverBase Inputs
            const double[::1] radius_array,
            const double[::1] density_array,
            const double[::1] gravity_array,
            const double complex[::1] shear_modulus_array,
            const double[::1] bulk_modulus_array,
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
            bool_cpp_t auto_solve = True
            ):

        # Loop variables
        cdef Py_ssize_t i

        # Load in floats and ints
        self.n_radius = len(radius_array)

        self.frequency  = frequency
        self.degree_l   = degree_l
        self.G_to_use   = G_to_use
        self.grav_coeff = 4. * pi * self.G_to_use

        # Load array information into class
        cdef np.ndarray[np.float64_t, ndim=1, mode='c'] radius_arr, bulk_arr, density_arr, gravity_arr
        cdef np.ndarray[np.complex128_t, ndim=1, mode='c'] shear_arr
        radius_arr  = np.empty(self.n_radius, dtype=np.float64, order='C')
        bulk_arr    = np.empty(self.n_radius, dtype=np.float64, order='C')
        density_arr = np.empty(self.n_radius, dtype=np.float64, order='C')
        gravity_arr = np.empty(self.n_radius, dtype=np.float64, order='C')
        shear_arr   = np.empty(self.n_radius, dtype=np.complex128, order='C')
        self.radius_view        = radius_arr
        self.shear_modulus_view = shear_arr
        self.bulk_modulus_view  = bulk_arr
        self.density_view       = density_arr
        self.gravity_view       = gravity_arr

        for i in range(self.n_radius):
            self.radius_view[i]        = radius_array[i]
            self.shear_modulus_view[i] = shear_modulus_array[i]
            self.bulk_modulus_view[i]  = bulk_modulus_array[i]
            self.density_view[i]       = density_array[i]
            self.gravity_view[i]       = gravity_array[i]

        # Initialize state variables
        self.shear_modulus = 0. + 0.j
        self.bulk_modulus  = 0.
        self.density = 0.
        self.gravity = 0.

        # Setup regular CySolver
        # Make sure to pass the radius array as the t_eval
        CySolver.__init__(
            self, t_span, y0, args, rtol, atol, max_step, first_step, rk_method,
            self.radius_view, capture_extra, num_extra, interpolate_extra, expected_size)

    cdef void update_constants(self) noexcept nogil:

        # Define degree l constants.
        cdef double degree_l_flt

        degree_l_flt = <double>self.degree_l
        self.lp1  = degree_l_flt + 1.
        self.lm1  = degree_l_flt - 1.
        self.llp1 = degree_l_flt * self.lp1

    cdef void update_interp(
            self,
            double radius,
            bool_cpp_t update_bulk = True,
            bool_cpp_t update_shear = True
            ) noexcept nogil:

        # Set state variables based on an interpolation using the provided radius.
        cdef (double, Py_ssize_t) interp_out
        cdef Py_ssize_t index_j

        # The first interpolation will be the slowest as it must find the closest index.
        # We will use this index in the other interpolations.
        interp_out = interpj(radius, self.radius_view, self.density_view)
        self.density = interp_out[0]
        index_j      = interp_out[1]

        self.gravity = interp(radius, self.radius_view, self.gravity_view, provided_j=index_j)
        if update_bulk:
            self.bulk_modulus = interp(radius, self.radius_view, self.bulk_modulus_view, provided_j=index_j)
        if update_shear:
            self.shear_modulus = interp_complex(radius, self.radius_view, self.shear_modulus_view, provided_j=index_j)
