cdef void cf_collapse_layer_solution(
    double complex* solution_ptr,
    double complex* constant_vector_ptr,
    double complex** storage_by_solution,
    double* layer_radius_ptr,
    double* layer_density_ptr,
    double* layer_gravity_ptr,
    double frequency_to_use,
    size_t layer_start_index,
    size_t num_layer_slices,
    unsigned char num_sols,
    unsigned char max_num_y,
    unsigned char num_ys,
    unsigned char num_output_ys,
    unsigned char ytype_i,
    int layer_type,
    bint layer_is_static,
    bint layer_is_incomp
    ) noexcept nogil
