from libcpp cimport bool as bool_cpp_t

cdef class TidalPyBaseExtensionClass:

    cdef str name_prefix
    cdef public str class_name
    cdef bool_cpp_t debug_mode
