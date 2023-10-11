from libcpp cimport bool as bool_cpp_t


cdef size_t cf_find_num_solutions(
        bool_cpp_t is_solid,
        bool_cpp_t is_static,
        bool_cpp_t is_incompressible
        ) noexcept nogil