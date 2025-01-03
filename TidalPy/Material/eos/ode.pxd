from libcpp cimport bool as cpp_bool

from CyRK cimport PreEvalFunc

cdef struct EOSOutput:
    double density
    double complex bulk_modulus
    double complex shear_modulus

cdef struct EOS_ODEInput:
    double G_to_use
    double planet_radius
    char* eos_input_ptr
    cpp_bool final_solve
    cpp_bool update_bulk
    cpp_bool update_shear

cdef void eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept nogil

