cdef void cf_get_surface_bc(
    double* boundary_conditions_ptr,
    size_t* bc_model_ptr,
    size_t num_bcs,
    double radius_to_use,
    double bulk_density_to_use,
    double degree_l_dbl,
    ) noexcept nogil
