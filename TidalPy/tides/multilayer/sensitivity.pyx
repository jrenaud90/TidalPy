# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
""" Functionality to calculate the sensitivity parameters described in Tobie et al (2005).

References
----------
TB05  : Tobie et al. (2005; DOI: 10.1016/j.icarus.2005.04.006)
KV21  : Kervazo et al. (2021; DOI: 10.1051/0004-6361/202039433)
"""

from libcpp cimport bool as cpp_bool

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.exceptions import ArgumentException
from TidalPy.constants cimport d_NAN_DBL
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx, cf_cabs, cf_cabs2
from TidalPy.RadialSolver.constants cimport MAX_NUM_Y


cdef void cf_calc_sensitivity_to_shear(
        double* radial_sensitivity_to_shear_ptr,
        double complex* radial_solutions_ptr, 
        double* radius_array_ptr,
        double complex* shear_modulus_array_ptr,
        double complex* bulk_modulus_array_ptr,
        size_t total_slices,
        size_t num_ytypes,
        int degree_l
        ) noexcept nogil:
    """ Calculates the radial sensitvity to shear dissipation.

    References
    ----------
    TB05: Eq. 33
    KV21: Eq. C.2

    """

    # Optimizations
    cdef double degree_l_dbl = <double>degree_l
    cdef double llp1         = degree_l_dbl * (degree_l_dbl + 1.)
    cdef double ll2m1lp2     = degree_l_dbl * (degree_l_dbl * degree_l_dbl - 1.) * (degree_l_dbl + 2.)

    # Extract the required radial solution values, isolating the real and imaginary portions.
    cdef double complex y1, y2, y3, y4
    cdef double complex y1_lower, y1_upper, y1_conj_rp1, y1_conj_rm1, y1_gradient_conj
    cdef double complex y1y3_term
    cdef double complex shear, bulk
    cdef double y1y3_term_abs2, y3_abs2, y4_abs2
    cdef double r, r2, dr0, dr1, dra, drb, drc
    cdef Py_ssize_t r_i
    cdef Py_ssize_t sgn_total_slices = <Py_ssize_t>total_slices
    cdef size_t ytype_i
    cdef size_t y_index
    cdef size_t num_output_ys = MAX_NUM_Y * num_ytypes
    for ytype_i in range(num_ytypes):
        for r_i in range(sgn_total_slices):
            y_index = r_i * num_output_ys + ytype_i * MAX_NUM_Y

            y1 = radial_solutions_ptr[y_index + 0]
            y2 = radial_solutions_ptr[y_index + 1]
            y3 = radial_solutions_ptr[y_index + 2]
            y4 = radial_solutions_ptr[y_index + 3]

            # Get y1 below and above this slice, used to find the derivative
            if r_i == 0:
                # Lower is undefined at center
                y1_lower = cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)
                y1_upper = radial_solutions_ptr[(r_i+1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
            elif r_i == (total_slices - 1):
                # Upper is undefined at surface
                y1_lower = radial_solutions_ptr[(r_i-1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
                y1_upper = cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)
            else:
                y1_lower = radial_solutions_ptr[(r_i-1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
                y1_upper = radial_solutions_ptr[(r_i+1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]

            # Shear and bulk may be real or complex
            shear = shear_modulus_array_ptr[r_i]
            bulk  = bulk_modulus_array_ptr[r_i]

            # Combine radial solutions into complex numbers
            y1_conj = cf_build_dblcmplx(y1.real, -y1.imag)

            # Find the gradient of y1_conj.
            #   The radius step size may not be even throughout the planet if there are
            #   layers with more slices. So we will use a 2nd order approach based on np.gradient method.
            r = radius_array_ptr[r_i]
            if r_i == 0:
                # First edge point
                y1_conj_rp1 = cf_build_dblcmplx(y1_upper.real, -y1_upper.imag)  # conj(y1) at r + 1
                dr1 = radius_array_ptr[r_i + 1] - r
                y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr1
            elif r_i == total_slices - 1:
                # Last end point
                y1_conj_rm1 = cf_build_dblcmplx(y1_lower.real, -y1_lower.imag) # conj(y1) at r - 1
                dr0 = r - radius_array_ptr[r_i - 1]
                y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
            else:
                # Midpoints
                y1_conj_rp1 = cf_build_dblcmplx(y1_upper.real, -y1_upper.imag)  # conj(y1) at r + 1
                y1_conj_rm1 = cf_build_dblcmplx(y1_lower.real, -y1_lower.imag)  # conj(y1) at r - 1
                dr0 = r - radius_array_ptr[r_i - 1]  # np.diff(r)[:-1]
                dr1 = radius_array_ptr[r_i + 1] - r  # np.diff(r)[0:]

                if dr0 == 0.0 and dr1 != 0.0:
                    # Treat as first edge point
                    y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr1
                elif dr0 != 0.0 and dr1 == 0.0:
                    # Treat as last end point
                    y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
                else:
                    # True midpoint
                    dra = -dr1 / (dr0 * (dr0 + dr1))
                    drb = (dr1 - dr0) / (dr0 * dr1)
                    drc = dr0 / (dr1 * (dr0 + dr1))
                    
                    y1_gradient_conj = dra * y1_conj_rm1 + drb * y1_conj + drc * y1_conj_rp1

            # Optimizations
            r2 = r * r
            y1y3_term      = 2.0 * y1 - llp1 * y3
            y1y3_term_abs2 = cf_cabs2(y1y3_term)
            y4_abs2 = cf_cabs2(y4)
            y3_abs2 = cf_cabs2(y3)

            # Find shear sensitivity
            if r == 0.:
                # The shear sensitivity is not defined at r=0
                radial_sensitivity_to_shear_ptr[r_i] = d_NAN_DBL
            else:
                radial_sensitivity_to_shear_ptr[r_i] = \
                    ((4. / 3.) * r2 / (cf_cabs2(bulk + (4. / 3.) * shear)) *
                    cf_cabs2(y2 - ((bulk - (2. / 3.) * shear) / r) * y1y3_term)) + \
                    (-(4. / 3.) * r * (y1_gradient_conj * y1y3_term).real + (1. / 3.) * y1y3_term_abs2) + \
                    ((llp1 * r2 * y4_abs2 / cf_cabs2(shear)) + (ll2m1lp2 * y3_abs2))

cdef void cf_calc_sensitivity_to_bulk(
        double* radial_sensitivity_to_bulk_ptr,
        double complex* radial_solutions_ptr, 
        double* radius_array_ptr,
        double complex* shear_modulus_array_ptr,
        double complex* bulk_modulus_array_ptr,
        size_t total_slices,
        size_t num_ytypes,
        int degree_l
        ) noexcept nogil:
    """ Calculates the radial sensitvity to bulk dissipation.

    References
    ----------
    TB05: Eq. 33
    KV21: Eq. C.1

    """

    # Optimizations
    cdef double degree_l_dbl = <double>degree_l
    cdef double llp1         = degree_l_dbl * (degree_l_dbl + 1.)
    cdef double ll2m1lp2     = degree_l_dbl * (degree_l_dbl * degree_l_dbl - 1.) * (degree_l_dbl + 2.)

    # Extract the required radial solution values, isolating the real and imaginary portions.
    cdef double complex y1, y2, y3, y4
    cdef double complex y1_lower, y1_upper, y1_conj_rp1, y1_conj_rm1, y1_gradient_conj
    cdef double complex y1y3_term
    cdef double complex shear, bulk
    cdef double y1y3_term_abs2, y3_abs2, y4_abs2
    cdef double r, r2, dr0, dr1, dra, drb, drc
    cdef Py_ssize_t r_i
    cdef Py_ssize_t sgn_total_slices = <Py_ssize_t>total_slices
    cdef size_t ytype_i
    cdef size_t y_index
    cdef size_t out_index
    cdef size_t num_output_ys = MAX_NUM_Y * num_ytypes
    for ytype_i in range(num_ytypes):
        out_index = ytype_i * total_slices
        for r_i in range(sgn_total_slices):
            y_index = r_i * num_output_ys + ytype_i * MAX_NUM_Y
            
            y1 = radial_solutions_ptr[y_index + 0]
            y2 = radial_solutions_ptr[y_index + 1]
            y3 = radial_solutions_ptr[y_index + 2]
    
            # Shear and bulk may be real or complex
            shear = shear_modulus_array_ptr[r_i]
            bulk  = bulk_modulus_array_ptr[r_i]
    
            # Combine radial solutions into complex numbers
            y1_conj = cf_build_dblcmplx(y1.real, -y1.imag)
    
            # Find the gradient of y1_conj.
            #   The radius step size may not be even throughout the planet if there are
            #   layers with more slices. So we will use a 2nd order approach based on np.gradient method.
            r = radius_array_ptr[r_i]
            if r_i == 0:
                # First edge point
                y1_upper = radial_solutions_ptr[(r_i+1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
                y1_conj_rp1 = cf_build_dblcmplx(y1_upper.real, -y1_upper.imag)  # conj(y1) at r + 1
                dr1 = radius_array_ptr[r_i + 1] - r
                y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr1
            elif r_i == total_slices - 1:
                # Last end point
                y1_lower = radial_solutions_ptr[(r_i-1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
                y1_conj_rm1 = cf_build_dblcmplx(y1_lower.real, -y1_lower.imag)  # conj(y1) at r - 1
                dr0 = r - radius_array_ptr[r_i - 1]
                y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
            else:
                # Midpoints
                y1_upper = radial_solutions_ptr[(r_i+1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
                y1_lower = radial_solutions_ptr[(r_i-1) * num_output_ys + ytype_i * MAX_NUM_Y + 0]
                y1_conj_rp1 = cf_build_dblcmplx(y1_upper.real, -y1_upper.imag)  # conj(y1) at r + 1
                y1_conj_rm1 = cf_build_dblcmplx(y1_lower.real, -y1_lower.imag)  # conj(y1) at r - 1
                dr0 = r - radius_array_ptr[r_i - 1]  # np.diff(r)[:-1]
                dr1 = radius_array_ptr[r_i + 1] - r  # np.diff(r)[0:]

                if dr0 == 0.0 and dr1 != 0.0:
                    # Treat as first edge point
                    y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr1
                elif dr0 != 0.0 and dr1 == 0.0:
                    # Treat as last end point
                    y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
                else:
                    # True midpoint
                    dra = -dr1 / (dr0 * (dr0 + dr1))
                    drb = (dr1 - dr0) / (dr0 * dr1)
                    drc = dr0 / (dr1 * (dr0 + dr1))
                    
                    y1_gradient_conj = dra * y1_conj_rm1 + drb * y1_conj + drc * y1_conj_rp1
    
            # Optimizations
            r2 = r * r
            y1y3_term = 2. * y1 - llp1 * y3
            y1y3_term_abs2 = cf_cabs2(y1y3_term)
    
            # Find bulk sensitivity
            if r == 0.:
                # The bulk sensitivity is not defined at r=0
                radial_sensitivity_to_bulk_ptr[out_index + r_i] = d_NAN_DBL
            else:
                radial_sensitivity_to_bulk_ptr[out_index + r_i] = \
                    (r2 / cf_cabs2(bulk + (4. / 3.) * shear)) * cf_cabs2(y2 - ((bulk - (2. / 3.) * shear) / r) * y1y3_term) + \
                    (2. * r * (y1_gradient_conj * y1y3_term).real) + \
                    y1y3_term_abs2


def calc_sensitivity_to_bulk(
        double complex[:, :] radial_solutions, 
        double[::1] radius_array,
        double complex[::1] shear_modulus_array,
        double complex[::1] bulk_modulus_array,
        int degree_l = 2,
        cpp_bool perform_checks = True
        ):

    cdef size_t total_slices = radius_array.size
    cdef size_t num_ytypes   = int(radial_solutions.shape[0] / 6)
    if perform_checks:
        if num_ytypes < 1:
            raise ArgumentException(f"Unexpected number of ytypes encountered ({num_ytypes}).")
        if shear_modulus_array.size != total_slices:
            raise ArgumentException("Unexpected size encountered for `shear_modulus_array`.")
        if bulk_modulus_array.size != total_slices:
            raise ArgumentException("Unexpected size encountered for `bulk_modulus_array`.")
        if radial_solutions.shape[1] != total_slices:
            raise ArgumentException("Unexpected size encountered for `radial_solutions` (axis=1).")

    # Build output array
    cdef cnp.ndarray sensitivity_to_bulk_array
    cdef double[::1] sensitivity_to_bulk_array_view
    if num_ytypes == 1:
        sensitivity_to_bulk_array = np.empty(total_slices, dtype=np.float64, order='C')
        sensitivity_to_bulk_array_view = sensitivity_to_bulk_array
    else:
        sensitivity_to_bulk_array = np.empty((num_ytypes, total_slices), dtype=np.float64, order='C')
        sensitivity_to_bulk_array_view = sensitivity_to_bulk_array.flatten()

    # Call cythonized function
    cf_calc_sensitivity_to_bulk(
        &sensitivity_to_bulk_array_view[0],
        &radial_solutions[0, 0], 
        &radius_array[0],
        &shear_modulus_array[0],
        &bulk_modulus_array[0],
        total_slices,
        num_ytypes,
        degree_l)

    return sensitivity_to_bulk_array

def calc_sensitivity_to_shear(
        double complex[:, :] radial_solutions, 
        double[::1] radius_array,
        double complex[::1] shear_modulus_array,
        double complex[::1] bulk_modulus_array,
        int degree_l = 2,
        cpp_bool perform_checks = True
        ):

    cdef size_t total_slices = radius_array.size
    cdef size_t num_ytypes   = int(radial_solutions.shape[0] / 6)
    if perform_checks:
        if num_ytypes < 1:
            raise ArgumentException(f"Unexpected number of ytypes encountered ({num_ytypes}).")
        if shear_modulus_array.size != total_slices:
            raise ArgumentException("Unexpected size encountered for `shear_modulus_array`.")
        if bulk_modulus_array.size != total_slices:
            raise ArgumentException("Unexpected size encountered for `bulk_modulus_array`.")
        if radial_solutions.shape[1] != total_slices:
            raise ArgumentException("Unexpected size encountered for `radial_solutions` (axis=1).")

    # Build output array
    cdef cnp.ndarray sensitivity_to_shear_array
    cdef double[::1] sensitivity_to_shear_array_view
    if num_ytypes == 1:
        sensitivity_to_shear_array = np.empty(total_slices, dtype=np.float64, order='C')
        sensitivity_to_shear_array_view = sensitivity_to_shear_array
    else:
        sensitivity_to_shear_array = np.empty((num_ytypes, total_slices), dtype=np.float64, order='C')
        sensitivity_to_shear_array_view = sensitivity_to_shear_array.flatten()

    # Call cythonized function
    cf_calc_sensitivity_to_shear(
        &sensitivity_to_shear_array_view[0],
        &radial_solutions[0, 0], 
        &radius_array[0],
        &shear_modulus_array[0],
        &bulk_modulus_array[0],
        total_slices,
        num_ytypes,
        degree_l)

    return sensitivity_to_shear_array
