cdef double complex cf_z_calc(
    double complex x_squared,
    int degree_l
    ) noexcept nogil

cdef void cf_takeuchi_phi_psi(
    double complex z,
    int degree_l,
    double complex* phi_ptr,
    double complex* phi_lplus1_ptr,
    double complex* psi_ptr,
    ) noexcept nogil