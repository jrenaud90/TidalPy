# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport fabs, fmax, isnan


cdef cpp_bool cf_isclose(
        double a,
        double b,
        double rtol = 1.0e-9,
        double atol = 0.0
        ) noexcept nogil:

    if isnan(a):
        return False
    
    if isnan(b):
        return False
    
    # Check for pure equivalence
    if a == b:
        return True

    # Check for closeness
    cdef double lhs = fabs(a - b)
    cdef double rhs = fmax(rtol * fmax(fabs(a), fabs(b)), atol)

    return lhs <= rhs


def isclose(
    double a,
    double b,
    double rtol = 1.0e-9,
    double atol = 0.0):

    return cf_isclose(a, b, rtol, atol)
