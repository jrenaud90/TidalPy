# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

def find_love(
    double complex[::1] complex_love_numbers_view,
    double complex[::1] surface_solutions_view,
    double surface_gravity
    ):
    """
    Find the complex Love and Shida numbers given the surface radial solutions for a planet.

    Parameters
    ----------
    complex_love_numbers_view : double complex[::1], array, output
        Array to store complex Love numbers. There must be space for 3 double complex numbers.
    surface_solutions_view : double complex[::1], array, input
        Array of radial solutions (y_i) values at the surface of a planet.
    surface_gravity : double, input
        Acceleration due to gravity at the planet's surface [m s-2].
    """

    # Create pointers; the c++ functions only work with doubles so we need to cast them to double
    cdef double* complex_love_numbers_ptr = <double*>&complex_love_numbers_view[0]
    cdef double* surface_solutions_ptr    = <double*>&surface_solutions_view[0]

    return find_love_cf(
        complex_love_numbers_ptr,
        surface_solutions_ptr,
        surface_gravity
        )
