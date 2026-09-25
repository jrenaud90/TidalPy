
cdef extern from "nondimensional_.hpp" nogil:
    cdef cppclass c_NonDimensionalScales:
        double second2_conversion
        double second_conversion
        double length_conversion
        double length3_conversion
        double density_conversion
        double mass_conversion
        double pascal_conversion

        c_NonDimensionalScales() except +
        c_NonDimensionalScales(
            double frequency,
            double mean_radius,
            double bulk_density
        )


cdef void cf_build_nondimensional_scales(
    c_NonDimensionalScales* non_dim_scales_ptr,
    double frequency,
    double mean_radius,
    double bulk_density
    ) noexcept nogil
