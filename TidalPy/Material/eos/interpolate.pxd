from libcpp cimport bool as cpp_bool

from TidalPy.Material.eos.ode cimport EOS_ODEInput

cdef struct InterpolateEOSInput:
    size_t num_slices
    double* radius_array_ptr
    double* density_array_ptr
    double complex* bulk_modulus_array_ptr
    double complex* shear_modulus_array_ptr

cdef void preeval_interpolate(
        # Values that will be updated by the function
        void* preeval_output,
        # Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        const void* preeval_input
        ) noexcept nogil
