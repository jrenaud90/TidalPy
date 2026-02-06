from libcpp cimport bool as cpp_bool

cdef extern from "numerics_.hpp" nogil:
    cpp_bool c_isclose(
        double a,
        double b,
        double rtol,
        double atol)