# distutils: language = c++

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
            double mean_radius,
            double bulk_density
        )
