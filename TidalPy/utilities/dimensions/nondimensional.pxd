
cdef extern from "nondimensional_.hpp" nogil:
    cdef cppclass NonDimensionalScalesCC:
        double second2_conversion
        double second_conversion
        double length_conversion
        double length3_conversion
        double density_conversion
        double mass_conversion
        double pascal_conversion


cdef void cf_build_nondimensional_scales(
    NonDimensionalScalesCC* non_dim_scales_ptr,
    double frequency,
    double mean_radius,
    double bulk_density
    ) noexcept nogil
