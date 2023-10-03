# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
Modified from https://github.com/numpy/numpy/blob/main/numpy/core/src/npymath/npy_math_complex.c.src */
/*==========================================================
 * Helper functions
 *
 * These are necessary because we do not count on using a
 * C99 compiler.
 *=========================================================*/
"""

from libc.math cimport NAN, INFINITY, isfinite, isnan, copysign, sqrt, fabs, signbit
from libc.float cimport DBL_MAX

# We risk spurious overflow for components >= DBL_MAX / (1 + sqrt(2)).
cdef double SQRT2 = 1.414213562373095048801688724209698079
cdef double THRESH = DBL_MAX / (1.0 + SQRT2)


cdef double hypot(double x, double y) noexcept nogil:

    cdef double yx
    cdef double temp

    if not isfinite(x) or not isfinite(y):
        return INFINITY

    if isnan(x) or isnan(y):
        return NAN

    x = fabs(x)
    y = fabs(y)
    if x < y:
        temp = x
        x    = y
        y    = temp

    if x == 0.:
        return 0.
    else:
        yx = y / x
        return x * sqrt(1. + yx * yx)

cdef double complex csqrt(double complex z) noexcept nogil:
    cdef double complex result
    cdef double a, b
    cdef double t
    cdef int scale

    a = z.real
    b = z.imag

    # Handle special cases.
    if b == 0.0:
        if a == 0.0:
            return 0.0 + 0.0j
        elif a > 0.0:
            return sqrt(a) + 0.0j

    if not isfinite(b):
        return INFINITY + 1.0j * b

    if isnan(a):
        # raise invalid if b is not a NaN
        t = (b - b) / (b - b)
        # Return NaN + NaN i
        return a + 1.0j * t

    if not isfinite(a):
        # csqrt(inf + NaN i)  = inf +  NaN i
        # csqrt(inf + y i)    = inf +  0 i
        # csqrt(-inf + NaN i) = NaN +- inf i
        # csqrt(-inf + y i)   = 0   +  inf i

        if signbit(a):
            return fabs(b - b) + 1.0j * copysign(a, b)
        else:
            return a + 1.0j * copysign(b - b, b)

    # The remaining special case (b is NaN) is handled just fine by the normal code path below.
    # Scale to avoid overflow.
    if fabs(a) >= THRESH or fabs(b) >= THRESH:
        a *= 0.25
        b *= 0.25
        scale = 1
    else:
        scale = 0

    # Algorithm 312, CACM vol 10, Oct 1967.
    if a >= 0.0:
        t = sqrt((a + hypot(a, b)) * 0.5)
        result = t + 1.0j * b / (2.0 * t)
    else:
        t = sqrt((-a + hypot(a, b)) * 0.5)
        result = fabs(b) / (2 * t) + 1.0j * copysign(t, b)

    # Rescale.
    if scale:
        return result.real * 2.0 + 1.0j * result.imag
    else:
        return result
