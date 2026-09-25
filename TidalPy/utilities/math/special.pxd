cdef double cf_double_factorial(int n) noexcept nogil

# ======================================================================================================================
# XSF Headers
# ======================================================================================================================
cdef extern from "xsf/bessel.h" namespace "xsf" nogil:
    # Needed for bessel functions used in sph_bessel.h
    pass

cdef extern from "xsf/sph_bessel.h" namespace "xsf" nogil:
    # Spherical Bessel j (first kind)
    double sph_bessel_j(long n, double x)                  # Real
    double complex sph_bessel_j(long n, double complex z)  # Complex
    double sph_bessel_j_jac(long n, double z)              # Derivative

    # Spherical Bessel y (second kind)
    double sph_bessel_y(long n, double x)                  # Real
    double complex sph_bessel_y(long n, double complex z)  # Complex
    double sph_bessel_y_jac(long n, double x)              # Derivative

    # Modified spherical Bessel i (first kind)
    double sph_bessel_i(long n, double x)                  # Real
    double complex sph_bessel_i(long n, double complex z)  # Complex
    double sph_bessel_i_jac(long n, double z)              # Derivative

    # Modified spherical Bessel k (second kind)
    double sph_bessel_k(long n, double z)                  # Real
    double complex sph_bessel_k(long n, double complex z)  # Complex
    double sph_bessel_k_jac(long n, double x)              # Derivative
