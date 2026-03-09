from libcpp.complex cimport complex as cpp_complex


cdef extern from "solid_matrix_.hpp" nogil:
    void c_fundamental_matrix(
        size_t first_slice_index,
        size_t num_radial_slices,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        cpp_complex[double]* complex_shear_array_ptr,
        cpp_complex[double]* fundamental_mtx_ptr,
        cpp_complex[double]* inverse_fundamental_mtx_ptr,
        cpp_complex[double]* derivative_mtx_ptr,
        int degree_l,
        double G_to_use
        ) noexcept
