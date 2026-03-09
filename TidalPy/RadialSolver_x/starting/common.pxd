from libcpp.complex cimport complex as cpp_complex


cdef extern from "common_.hpp" nogil:

    cdef cpp_complex[double] c_spherical_jn(
        int n,
        cpp_complex[double] z) noexcept nogil

    cdef cpp_complex[double] c_z_calc(
        cpp_complex[double] x_squared,
        int degree_l) noexcept nogil

    cdef void c_takeuchi_phi_psi(
        cpp_complex[double] z2,
        int degree_l,
        cpp_complex[double]* phi_ptr,
        cpp_complex[double]* phi_lplus1_ptr,
        cpp_complex[double]* psi_ptr) noexcept nogil
