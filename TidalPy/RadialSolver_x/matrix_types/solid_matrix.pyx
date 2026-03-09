# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.RadialSolver_x.matrix_types.solid_matrix cimport c_fundamental_matrix


def fundamental_matrix(
        double[::1] radius_array_view,
        double[::1] density_array_view,
        double[::1] gravity_array_view,
        double complex[::1] complex_shear_array_view,
        int degree_l = 2,
        double G_to_use = 6.67430e-11,
        cpp_bool perform_checks = True
        ):
    """ Construct the fundamental matrix and its inverse using harmonic degree l.

    See Eq. 2.42 of SVC16

    Assumptions
    -----------
    - These matrices assume an incompressible body.

    Parameters
    ----------
    radius_array_view : double[::1]
        Array of radius values [m].
    density_array_view : double[::1]
        Array of density at each radius [kg m-3].
    gravity_array_view : double[::1]
        Array of acceleration due to gravity at each radius [m s-2].
    complex_shear_array_view : double complex[::1]
        Array of complex shear modulus at each radius [Pa].
    degree_l : int, default=2
        Harmonic degree.
    G_to_use : double, default=6.67430e-11
        Gravitational constant.
    perform_checks : bool, default=True
        If True, checks will be performed on input arguments.
    """

    cdef size_t num_radial_slices = radius_array_view.size

    if perform_checks:
        if <size_t>density_array_view.size != num_radial_slices:
            raise ValueError('Unexpected size encountered for density array.')
        if <size_t>gravity_array_view.size != num_radial_slices:
            raise ValueError('Unexpected size encountered for gravity array.')
        if <size_t>complex_shear_array_view.size != num_radial_slices:
            raise ValueError('Unexpected size encountered for complex shear array.')

    # Build output arrays
    cdef cnp.ndarray[cnp.complex128_t, ndim=3] fundamental_mtx_arr         = np.empty((num_radial_slices, 6, 6), dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=3] inverse_fundamental_mtx_arr = np.empty((num_radial_slices, 6, 6), dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=3] derivative_mtx_arr          = np.empty((num_radial_slices, 6, 6), dtype=np.complex128, order='C')

    cdef double complex[:, :, ::1] fundamental_mtx_view         = fundamental_mtx_arr
    cdef double complex[:, :, ::1] inverse_fundamental_mtx_view = inverse_fundamental_mtx_arr
    cdef double complex[:, :, ::1] derivative_mtx_view          = derivative_mtx_arr

    c_fundamental_matrix(
        0,
        num_radial_slices,
        &radius_array_view[0],
        &density_array_view[0],
        &gravity_array_view[0],
        <cpp_complex[double]*>&complex_shear_array_view[0],
        <cpp_complex[double]*>&fundamental_mtx_view[0, 0, 0],
        <cpp_complex[double]*>&inverse_fundamental_mtx_view[0, 0, 0],
        <cpp_complex[double]*>&derivative_mtx_view[0, 0, 0],
        degree_l,
        G_to_use
    )

    return fundamental_mtx_arr, inverse_fundamental_mtx_arr, derivative_mtx_arr
