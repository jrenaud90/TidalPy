from libcpp.complex cimport complex as cpp_complex


cdef extern from "takeuchi_.hpp" nogil:

    cdef void c_takeuchi_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        cpp_complex[double] shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_takeuchi_solid_static_compressible(
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        cpp_complex[double] shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_takeuchi_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil
