cdef void cf_calc_sensitivity_to_shear(
        double* radial_sensitivity_to_shear_ptr,
        double complex* radial_solutions_ptr, 
        double* radius_array_ptr,
        double complex* shear_modulus_array_ptr,
        double complex* bulk_modulus_array_ptr,
        size_t total_slices,
        size_t num_ytypes,
        int degree_l
        ) noexcept nogil

cdef void cf_calc_sensitivity_to_bulk(
    double* radial_sensitivity_to_bulk_ptr,
    double complex* radial_solutions_ptr, 
    double* radius_array_ptr,
    double complex* shear_modulus_array_ptr,
    double complex* bulk_modulus_array_ptr,
    size_t total_slices,
    size_t num_ytypes,
    int degree_l
    ) noexcept nogil
