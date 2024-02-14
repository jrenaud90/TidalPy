# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.utilities.math.complex cimport cf_build_dblcmplx

cdef void find_love_cf(
        double complex* complex_love_numbers_ptr,
        double complex* surface_solutions_ptr,
        double surface_gravity) noexcept nogil:

    # Extract the required radial solution values, isolating the real and imaginary portions.
    cdef double complex y1 = surface_solutions_ptr[0]
    cdef double complex y3 = surface_solutions_ptr[2]
    cdef double complex y5 = surface_solutions_ptr[4]
    cdef double complex c_surface_gravity = cf_build_dblcmplx(surface_gravity, 0.0)

    # Calculate Love and Shida numbers
    # Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), opposite convention of Tobie et al. (2005, eqs. 9 & 36)
    # And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9.
     
    # Okay, some clarification on this. It looks like VS04 that HH14 is based on used a different convention for y5,
    #      Tobie05's y5 = -y5 of SV; we follow that format here.
    
    # Love k
    complex_love_numbers_ptr[0] = cf_build_dblcmplx(y5.real - 1.0, y5.imag)
    # Love h
    complex_love_numbers_ptr[1] = y1 * c_surface_gravity
    # Shida l
    complex_love_numbers_ptr[2] = y3 * c_surface_gravity


def find_love(
    double complex[:] complex_love_numbers_view,
    double complex[:] surface_solutions_view,
    double surface_gravity
    ):
    """
    Find the complex Love and Shida numbers given the surface radial solutions for a planet.

    Parameters
    ----------
    complex_love_numbers_view : double complex[:], array, output
        Array to store complex Love numbers. There must be space for 3 double complex numbers.
    surface_solutions_view : double complex[:], array, input
        Array of radial solutions (y_i) values at the surface of a planet.
    surface_gravity : double, input
        Acceleration due to gravity at the planet's surface [m s-2].
    """

    return find_love_cf(
        &complex_love_numbers_view[0],
        &surface_solutions_view[0],
        surface_gravity
        )
