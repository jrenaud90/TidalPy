from libcpp.complex cimport complex as cpp_complex


cdef extern from "kamata_.hpp" nogil:

    cdef void c_kamata_solid_dynamic_compressible(
        const double frequency,
        const double radius,
        const double density,
        const cpp_complex[double]& bulk_modulus,
        const cpp_complex[double]& shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_solid_static_compressible(
        const double radius,
        const double density,
        const cpp_complex[double]& bulk_modulus,
        const cpp_complex[double]& shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_solid_dynamic_incompressible(
        const double frequency,
        const double radius,
        const double density,
        const cpp_complex[double]& shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_liquid_dynamic_compressible(
        const double frequency,
        const double radius,
        const double density,
        const cpp_complex[double]& bulk_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef void c_kamata_liquid_dynamic_incompressible(
        const double frequency,
        const double radius,
        const double density,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil
