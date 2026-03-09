from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex


cdef extern from "collapse_.hpp" nogil:

    cdef void c_collapse_layer_solution(
        cpp_complex[double]* solution_ptr,
        cpp_complex[double]* constant_vector_ptr,
        cpp_complex[double]** storage_by_solution,
        double* layer_radius_ptr,
        double* layer_density_ptr,
        double* layer_gravity_ptr,
        double frequency_to_use,
        size_t layer_start_index,
        size_t num_layer_slices,
        size_t num_sols,
        size_t max_num_y,
        size_t num_ys,
        size_t num_output_ys,
        size_t ytype_i,
        int layer_type,
        cpp_bool layer_is_static,
        cpp_bool layer_is_incomp) noexcept nogil
