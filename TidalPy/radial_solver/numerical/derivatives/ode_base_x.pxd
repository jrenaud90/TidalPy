from libcpp cimport bool as bool_cpp_t

from CyRK.cy.cysolver cimport CySolver


cdef class RadialSolverBase(CySolver):

    # Global constants
    cdef double frequency
    cdef unsigned int degree_l
    cdef double degree_l_flt, lp1, lm1, llp1
    cdef double G_to_use, grav_coeff

    # Number of radial slices
    cdef Py_ssize_t num_slices

    # Radial pointer arrays
    cdef double* radius_array_ptr
    cdef double* density_array_ptr
    cdef double* gravity_array_ptr
    cdef double* bulk_modulus_array_ptr
    cdef double complex* shear_modulus_array_ptr

    # State properties at current radius
    cdef double complex shear_modulus
    cdef double bulk_modulus
    cdef double density
    cdef double gravity

    # Class methods
    cdef void install_pointers(
        self,

        # RadialSolverBase pointers
        Py_ssize_t num_slices,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_modulus_array_ptr,
        double complex* shear_modulus_array_ptr,
        double* atols,
        double* rtols,

        # Additional optional arguments for RadialSolver class
        bool_cpp_t limit_solution_to_radius = *,
        bool_cpp_t call_first_reset = *,
        bool_cpp_t auto_solve = *
        )

    cdef void update_interp(
            self,
            bool_cpp_t update_bulk,
            bool_cpp_t update_shear
            ) noexcept nogil
