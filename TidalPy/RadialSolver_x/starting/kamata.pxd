from libcpp.complex cimport complex as cpp_complex


cdef extern from "kamata_.hpp" nogil:

    cdef void c_kamata_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        cpp_complex[double] shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_solid_static_compressible(
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        cpp_complex[double] shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_solid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        cpp_complex[double] shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        cpp_complex[double] bulk_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_liquid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        int degree_l,
        double G_to_use,
        size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil
