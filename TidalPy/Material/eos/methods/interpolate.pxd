from libcpp cimport bool as cpp_bool

from TidalPy.Material.eos.ode cimport EOS_ODEInput

cdef int EOS_INTERPOLATE_METHOD_INT = 0

cdef struct InterpolateEOSInput:
    size_t num_slices
    double* radius_array_ptr
    double* density_array_ptr
    double complex* bulk_modulus_array_ptr
    double complex* shear_modulus_array_ptr

cdef void preeval_interpolate(
        # Values that will be updated by the function
        char* preeval_output,
        # Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept nogil
