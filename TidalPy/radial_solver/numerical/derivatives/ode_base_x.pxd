from libcpp cimport bool as bool_cpp_t

from CyRK.cy.cysolver cimport CySolver


cdef class RadialSolverBase(CySolver):

    # Class Attributes
    cdef unsigned int n_radius

    # Global constants
    cdef double frequency
    cdef unsigned int degree_l
    cdef double degree_l_flt, lp1, lm1, llp1
    cdef double G_to_use, grav_coeff

    # Radial arrays
    cdef double[::1] radius_view,
    cdef double complex[::1] shear_modulus_view,
    cdef double[::1] bulk_modulus_view,
    cdef double[::1] density_view,
    cdef double[::1] gravity_view,

    # State properties at current radius
    cdef double complex shear_modulus
    cdef double bulk_modulus
    cdef double density
    cdef double gravity

    # Class methods
    cdef void update_interp(self, bool_cpp_t update_bulk, bool_cpp_t update_shear) noexcept nogil
