from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


cdef extern from "boundaries_.hpp" nogil:

    cdef void c_apply_surface_bc(
        cpp_complex[double]* constant_vector_ptr,
        int* bc_solution_info_ptr,
        double* bc_pointer,
        cpp_complex[double]* uppermost_y_per_solution_ptr,
        double surface_gravity,
        double G_to_use,
        size_t num_sols,
        size_t max_num_y,
        size_t ytype_i,
        int layer_type,
        cpp_bool layer_is_static,
        cpp_bool layer_is_incomp) noexcept nogil
