cdef void cf_fundamental_matrix(
    Py_ssize_t num_radial_slices,
    double* radius_array_ptr,
    double* density_array_ptr,
    double* gravity_array_ptr,
    double complex* complex_shear_array_ptr,
    double complex* fundamental_mtx_ptr,
    double complex* inverse_fundamental_mtx_ptr,
    double complex* derivative_mtx_ptr,
    unsigned char degree_l = *,
    double G_to_use = *
    ) noexcept nogil
