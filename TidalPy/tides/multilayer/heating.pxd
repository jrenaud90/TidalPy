
cdef void cf_calc_radial_volumetric_tidal_heating(
    double* volumetric_tidal_heating_arr_ptr,
    size_t total_slices,
    double eccentricity,
    double orbital_frequency,
    double semi_major_axis,
    double host_mass,
    double* radius_arr_ptr,
    double* radial_sensitivity_to_shear_arr_ptr,
    double complex* complex_shear_modulus_arr_ptr,
    double* radial_sensitivity_to_bulk_arr_ptr,
    double complex* complex_bulk_modulus_arr_ptr,
    int degree_l,
    double G_to_use
    ) noexcept nogil
