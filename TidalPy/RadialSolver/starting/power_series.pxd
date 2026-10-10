from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


cdef extern from "power_series_.hpp" nogil:

    cdef cpp_bool c_power_series_solid(
        const double frequency,
        const double radius,
        const double density,
        const cpp_complex[double]& bulk_modulus,
        const cpp_complex[double]& shear_modulus,
        const cpp_bool is_incompressible,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil

    cdef cpp_bool c_power_series_liquid_dynamic(
        const double frequency,
        const double radius,
        const double density,
        const cpp_complex[double]& bulk_modulus,
        const cpp_bool is_incompressible,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        cpp_complex[double]* starting_conditions_ptr) noexcept nogil
