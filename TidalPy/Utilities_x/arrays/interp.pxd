# distutils: language = c++
"""Cython declarations for TidalPy's 1-D linear interpolation utilities (``interp_.hpp``)."""

from libcpp.complex cimport complex as cpp_complex


cdef extern from "interp_.hpp" namespace "tidalpy" nogil:
    double c_interp(
        double desired_x,
        const double* x_domain,
        const double* dependent_values,
        size_t len_x,
        size_t guess)

    cpp_complex[double] c_interp_complex(
        double desired_x,
        const double* x_domain,
        const cpp_complex[double]* dependent_values,
        size_t len_x,
        size_t guess)
