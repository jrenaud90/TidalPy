from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

from CyRK cimport PreEvalFunc


cdef extern from "ode_.hpp" nogil:

    const size_t C_EOS_Y_VALUES
    const size_t C_EOS_EXTRA_VALUES
    const size_t C_EOS_DY_VALUES

    cdef struct c_EOSOutput:
        double density
        cpp_complex[double] bulk_modulus
        cpp_complex[double] shear_modulus

    cdef struct c_EOS_ODEInput:
        double G_to_use
        double planet_radius
        char*  eos_input_ptr
        cpp_bool final_solve
        cpp_bool update_bulk
        cpp_bool update_shear

    void c_eos_diffeq(
            double* dy_ptr,
            double radius,
            double* y_ptr,
            char* input_args,
            PreEvalFunc eos_function) noexcept
