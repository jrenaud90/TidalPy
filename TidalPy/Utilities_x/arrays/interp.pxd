# distutils: language = c++
"""
interp.pxd
Cython declarations for TidalPy's 1-D linear interpolation utilities.

Other extensions can cimport these to interpolate at C speed::

    from TidalPy.Utilities_x.arrays.interp cimport c_interp, c_interp_complex
"""

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
