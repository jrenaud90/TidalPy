# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp.complex cimport complex as cpp_complex
from cpython.complex cimport PyComplex_FromDoubles


def z_calc(double complex x_squared, int degree_l):
    """
    Calculate the z function using spherical Bessel function.

    See Eq. B14 of KMN15, Eqs. 96-97 of TS72.

    Parameters
    ----------
    x_squared : complex
        Expression passed to the Bessel function.
    degree_l : int
        Tidal harmonic order.

    Returns
    -------
    z : complex
    """
    cdef cpp_complex[double] x_sq = cpp_complex[double](x_squared.real, x_squared.imag)
    cdef cpp_complex[double] result = c_z_calc(x_sq, degree_l)
    return PyComplex_FromDoubles(result.real(), result.imag())


def takeuchi_phi_psi(double complex z2, int degree_l):
    """
    Calculate phi, phi_{l+1}, and psi functions for Takeuchi starting conditions.

    See TS72 Eq. 103.

    Parameters
    ----------
    z2 : complex
        z^2 argument.
    degree_l : int
        Tidal harmonic order.

    Returns
    -------
    phi : complex
    phi_lplus1 : complex
    psi : complex
    """
    cdef cpp_complex[double] z2_cpp = cpp_complex[double](z2.real, z2.imag)
    cdef cpp_complex[double] phi, phi_lp1, psi
    c_takeuchi_phi_psi(z2_cpp, degree_l, &phi, &phi_lp1, &psi)
    return (
        PyComplex_FromDoubles(phi.real(), phi.imag()),
        PyComplex_FromDoubles(phi_lp1.real(), phi_lp1.imag()),
        PyComplex_FromDoubles(psi.real(), psi.imag()),
    )
