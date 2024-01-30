from TidalPy.utilities.types_x cimport double_numeric

cdef void cf_non_dimensionalize_physicals(
        size_t num_radius,
        double frequency,
        double mean_radius,
        double bulk_density,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_array_ptr,
        double_numeric* shear_array_ptr,
        double* frequency_to_use,
        double* G_to_use
        ) noexcept nogil

cdef void cf_redimensionalize_physicals(
        size_t num_radius,
        double frequency,
        double mean_radius,
        double bulk_density,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_array_ptr,
        double_numeric* shear_array_ptr,
        double* frequency_to_use,
        double* G_to_use
        ) noexcept nogil

cdef void cf_redimensionalize_radial_functions(
    double complex* radial_function_ptr,
    double mean_radius,
    double bulk_density,
    size_t num_slices,
    size_t num_solutions = *) noexcept nogil
