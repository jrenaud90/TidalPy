from libcpp cimport bool as cpp_bool

cdef cpp_bool cf_isclose(
    double a,
    double b,
    double rtol = *,
    double atol = *
    ) noexcept nogil