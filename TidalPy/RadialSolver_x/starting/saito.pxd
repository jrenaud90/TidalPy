from libcpp.complex cimport complex as cpp_complex


cdef extern from "saito_.hpp" nogil:

    cdef void c_saito_liquid_static_incompressible(
        const double radius,
        const int degree_l,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil
