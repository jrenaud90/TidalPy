from libcpp cimport bool as bool_cpp_t

from CyRK.cy.cysolver cimport CySolver


cdef class RadialSolverBase(CySolver):

    # Global constants
    cdef double frequency
    cdef unsigned int degree_l
    cdef double degree_l_flt, lp1, lm1, llp1
    cdef double G_to_use, grav_coeff

    # Radial arrays
    cdef Py_ssize_t n_radius
    cdef double* radius_array_ptr
    cdef double complex* shear_modulus_array_ptr
    cdef double* bulk_modulus_array_ptr
    cdef double* density_array_ptr
    cdef double* gravity_array_ptr

    # State properties at current radius
    cdef double complex shear_modulus
    cdef double bulk_modulus
    cdef double density
    cdef double gravity

    # Class methods
    cdef void update_interp(
            self,
            bool_cpp_t update_bulk,
            bool_cpp_t update_shear
            ) noexcept nogil
