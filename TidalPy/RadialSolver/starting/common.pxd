cdef double complex cf_z_calc(
    double complex x_squared,
    unsigned char degree_l
    ) noexcept nogil

cdef void cf_takeuchi_phi_psi(
    double complex z,
    unsigned char degree_l,
    double complex* phi_ptr,
    double complex* phi_lplus1_ptr,
    double complex* psi_ptr,
    ) noexcept nogil