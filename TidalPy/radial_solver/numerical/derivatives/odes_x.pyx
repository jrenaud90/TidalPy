# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as np
np.import_array()

from libcpp cimport bool as bool_cpp_t

from CyRK.cy.cysolver cimport CySolver, MAX_STEP
from CyRK.array.interp cimport interpj, interp, interp_complex

from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_dynamic_incomp_x cimport dy_solid_dynamic_incompressible_x, dy_liquid_dynamic_incompressible_x
from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_dynamic_x cimport dy_solid_dynamic_compressible_x, dy_liquid_dynamic_compressible_x
from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_static_incomp_x cimport dy_solid_static_incompressible_x, dy_liquid_static_incompressible_x
from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_static_x cimport dy_solid_static_compressible_x, dy_liquid_static_compressible_x


cdef class BaseODE(CySolver):

    cdef unsigned int n_radius

    cdef double frequency
    cdef unsigned int degree_l
    cdef double G_to_use

    cdef double[:] radius_view,
    cdef double complex[:] shear_modulus_view,
    cdef double[:] bulk_modulus_view,
    cdef double[:] density_view,
    cdef double[:] gravity_view,

    cdef double complex shear
    cdef double bulk
    cdef double density
    cdef double gravity


    def __init__(
            self,
            const double[:] radius_array,
            const double complex[:] shear_modulus_array,
            const double[:] bulk_modulus_array,
            const double[:] density_array,
            const double[:] gravity_array,
            double frequency,
            unsigned int degree_l,
            double G_to_use,

            # Regular CySolver Inputs
            (double, double) t_span,
            const double[:] y0,
            tuple args = None,
            double rtol = 1.e-7,
            double atol = 1.e-8,
            double max_step = MAX_STEP,
            double first_step = 0.,
            unsigned char rk_method = 1,
            const double[:] t_eval = None,
            bool_cpp_t capture_extra = False,
            unsigned short num_extra = 0,
            bool_cpp_t interpolate_extra = False,
            unsigned int expected_size = 0
            ):

        # Loop variables
        cdef Py_ssize_t i

        # Load in floats and ints
        self.n_radius = len(radius_array)

        self.frequency = frequency
        self.degree_l  = degree_l
        self.G_to_use  = G_to_use

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
        self.shear   = 0. + 0.j
        self.bulk    = 0.
        self.density = 0.
        self.gravity = 0.

        # Setup regular CySolver
        # Make sure to pass the radius array as the t_eval
        CySolver.__init__(
            self, t_span, y0, args, rtol, atol, max_step, first_step, rk_method,
            self.radius_view, capture_extra, num_extra, interpolate_extra, expected_size)

    cdef void update_interp(
            self,
            double radius,
            bool_cpp_t update_bulk = True,
            bool_cpp_t update_shear = True,
            ) nogil:

        # Set state variables based on an interpolation using the provided radius.
        cdef (double, Py_ssize_t) first_out
        cdef Py_ssize_t index_j

        # The first interpolation will be the slowest as it must find the closest index.
        # We will use this index in the other interpolations.
        first_out = interpj(radius, self.radius_view, self.density_view)
        self.density = first_out[0]
        index_j      = first_out[1]

        self.gravity = interp(radius, self.radius_view, self.gravity_view, provided_j=index_j)
        if update_bulk:
            self.bulk = interp(radius, self.radius_view, self.bulk_modulus_view, provided_j=index_j)
        if update_shear:
            self.shear = interp_complex(radius, self.radius_view, self.shear_modulus_view, provided_j=index_j)


cdef class SolidDynamicCompressible(BaseODE):
    cdef void diffeq(self):

        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new)

        # Call the correct differential equation
        dy_solid_dynamic_compressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.shear, self.bulk, self.density, self.gravity, self.frequency,self.degree_l, self.G_to_use
            )


cdef class SolidDynamicIncompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new, update_bulk=False)

        # Call the correct differential equation
        dy_solid_dynamic_incompressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.shear, self.density, self.gravity, self.frequency, self.degree_l, self.G_to_use
            )


cdef class SolidStaticCompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new)

        # Call the correct differential equation
        dy_solid_static_compressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.shear, self.bulk, self.density, self.gravity, self.degree_l, self.G_to_use
            )


cdef class SolidStaticIncompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new, update_bulk=False)

        # Call the correct differential equation
        dy_solid_static_incompressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.shear, self.density, self.gravity, self.degree_l, self.G_to_use
            )


cdef class LiquidDynamicCompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new, update_bulk=True, update_shear=False)

        # Call the correct differential equation
        dy_liquid_dynamic_compressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.bulk, self.density, self.gravity, self.frequency, self.degree_l, self.G_to_use
            )


cdef class LiquidDynamicIncompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new, update_bulk=False, update_shear=False)

        # Call the correct differential equation
        dy_liquid_dynamic_incompressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.density, self.gravity, self.frequency, self.degree_l, self.G_to_use
            )


cdef class LiquidStaticCompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new, update_bulk=False, update_shear=False)

        # Call the correct differential equation
        dy_liquid_static_compressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.density, self.gravity, self.degree_l, self.G_to_use
            )


cdef class LiquidStaticIncompressible(BaseODE):
    cdef void diffeq(self):
        # Note that `t_new` is the current "radius" not time.

        # Update interpolation
        self.update_interp(self.t_new, update_bulk=False, update_shear=False)

        # Call the correct differential equation
        dy_liquid_static_incompressible_x(
            self.t_new, self.y_new_view, self.dy_new_view,
            self.density, self.gravity, self.degree_l, self.G_to_use
            )