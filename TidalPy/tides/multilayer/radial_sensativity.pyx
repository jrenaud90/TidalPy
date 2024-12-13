# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
""" Functionality to calculate the sensitivity parameters described in Tobie et al (2005).

References
----------
TB05  : Tobie et al. (2005; DOI: 10.1016/j.icarus.2005.04.006)
KV21  : Kervazo et al. (2021; DOI: 10.1051/0004-6361/202039433)
"""

from libc.math cimport abs

from cython.parallel cimport prange

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.constants cimport d_NAN_DBL
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx, cf_cabs, cf_cabs2
from TidalPy.RadialSolver.constants cimport MAX_NUM_Y


cdef cf_sensitivity_to_shear(
        double complex* radial_sensitivity_to_shear_ptr,
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
    cdef double y1y3_term_abs2, y3_abs2, y4_abs2
    cdef double dr0, dr1, dra, drb, drc, r_2
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
            y2 = y2.real + 1.0j * y2.imag
            y1_conj = y1.real - 1.0j * y1.imag

            # Find the gradient of y1_conj.
            #   The radius step size may not be even throughout the planet if there are
            #   layers with more slices. So we will use a 2nd order approach based on np.gradient method.
            r = radius_array_ptr[r_i]
            if r_i == 0:
                # First edge point
                y1_conj_rp1 = y1_upper.real - 1.0j * y1_upper.imag  # conj(y1) at r + 1
                dr0 = radius_array_ptr[r_i + 1] - r
                y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr0
            elif r_i == total_slices - 1:
                # Last end point
                y1_conj_rm1 = y1_lower.real - 1.0j * y1_lower.imag # conj(y1) at r - 1
                dr0 = r - radius_array_ptr[r_i - 1]
                y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
            else:
                # Midpoints
                y1_conj_rp1 = y1_upper.real - 1.0j * y1_upper.imag  # conj(y1) at r + 1
                y1_conj_rm1 = y1_lower.real - 1.0j * y1_lower.imag  # conj(y1) at r - 1
                dr0 = r - radius_array_ptr[r_i - 1]  # np.diff(r)[:-1]
                dr1 = radius_array_ptr[r_i + 1] - r  # np.diff(r)[0:]

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
                radial_sensitivity_to_shear_ptr[r_i] = cf_build_dblcmplx(d_NAN_DBL, d_NAN_DBL)
            else:
                radial_sensitivity_to_shear_ptr[r_i] = \
                    ((4. / 3.) * r2 / (cf_cabs2(bulk + (4. / 3.) * shear)) *
                    cf_cabs2(y2 - ((bulk - (2. / 3.) * shear) / r) * y1y3_term)) + \
                    (-(4. / 3.) * r * (y1_gradient_conj * y1y3_term).real + (1. / 3.) * y1y3_term_abs2) + \
                    ((llp1 * r2 * y4_abs2 / (cf_cabs(shear)**2)) + (ll2m1lp2 * y3_abs2))

    return radial_sensitivity_to_shear_ptr

@njit(cacheable=True)
def sensitivity_to_bulk(
        radial_solutions_ptr: np.ndarray, radius_array_ptr: np.ndarray,
        shear_modulus_array_ptr: np.ndarray, bulk_modulus_array_ptr: np.ndarray, order_l: int = 2
        ):
    """ Calculates the radial sensitvity to bulk dissipation.

    References
    ----------
    TB05: Eq. 33
    KV21: Eq. C.1

    Parameters
    ----------
    radial_solutions_ptr : np.ndarray
        Viscoelastic-Gravitational radial solutions found through either the matrix propagation or
        numerical integration techniques.
        These should follow the TS72 order convention and dimensions. With the exception that there are double
        the number of values to account for the imaginary portions. E.g., y1.real, y1.imag, ... (Matrix: 12 x N)
    radius_array_ptr : np.ndarray
        Radius at the top of each radial slice throughout the world (length N) [m]
    shear_modulus_array_ptr : np.ndarray
        The shear modulus defined at each radius slice (length N) [Pa].
        Can be real or complex valued.
    bulk_modulus_array_ptr : np.ndarray
        The bulk modulus defined at each radius slice (length N) [Pa].
        Can be real or complex valued.
        Note from ID:
            "K=200.0e9; // Bulk modulus (g cm-2 s-2), arbitrarily higher than
            K_rock (39-133 GPa) and K_ice (10.7 GPa) for consistency with incompressible prop mtx."
    order_l : int = 2
        Tidal Harmonic Degree

    Returns
    -------
    radial_sensitivity_to_bulk : np.ndarray
        Radial sensitivity to shear stress as defined in TB05, Eq. 33. This is a real-valued float-array (length N)

    """

    # Basic information
    total_slices = radius_array_ptr.shape[0]

    # Build arrays
    radial_sensitivity_to_bulk = np.empty(total_slices, dtype=np.float64)

    # Optimizations
    llp1     = order_l * (order_l + 1.)

    # Extract the required radial solution values, isolating the real and imaginary portions.
    for r_i in range(total_slices):
        y1 = radial_solutions_ptr[0, r_i]
        y2 = radial_solutions_ptr[1, r_i]
        y3 = radial_solutions_ptr[2, r_i]

        # Shear and bulk may be real or complex
        shear = shear_modulus_array_ptr[r_i]
        bulk  = bulk_modulus_array_ptr[r_i]

        # Combine radial solutions into complex numbers
        y2 = y2.real + 1.0j * y2.imag
        y1_conj = y1.real - 1.0j * y1.imag

        # Find the gradient of y1_conj.
        #   The radius step size may not be even throughout the planet if there are
        #   layers with more slices. So we will use a 2nd order approach based on np.gradient method.
        r = radius_array_ptr[r_i]
        if r_i == 0:
            # First edge point
            y1_upper.real = np.real(radial_solutions_ptr[0, r_i + 1])
            y1_upper.imag = np.imag(radial_solutions_ptr[0, r_i + 1])
            y1_conj_rp1 = y1_upper.real - 1.0j * y1_upper.imag  # conj(y1) at r + 1
            dr0 = radius_array_ptr[r_i + 1] - r
            y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr0
        elif r_i == total_slices - 1:
            # Last end point
            y1_lower.real = np.real(radial_solutions_ptr[0, r_i - 1])
            y1_lower.imag = np.imag(radial_solutions_ptr[0, r_i - 1])
            y1_conj_rm1 = y1_lower.real - 1.0j * y1_lower.imag  # conj(y1) at r - 1
            dr0 = r - radius_array_ptr[r_i - 1]
            y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
        else:
            # Midpoints
            y1_upper.real = np.real(radial_solutions_ptr[0, r_i + 1])
            y1_upper.imag = np.imag(radial_solutions_ptr[0, r_i + 1])
            y1_lower.real = np.real(radial_solutions_ptr[0, r_i - 1])
            y1_lower.imag = np.imag(radial_solutions_ptr[0, r_i - 1])
            y1_conj_rp1 = y1_upper.real - 1.0j * y1_upper.imag  # conj(y1) at r + 1
            y1_conj_rm1 = y1_lower.real - 1.0j * y1_lower.imag  # conj(y1) at r - 1
            dr0 = r - radius_array_ptr[r_i - 1]  # np.diff(r)[:-1]
            dr1 = radius_array_ptr[r_i + 1] - r  # np.diff(r)[0:]

            dra = -dr1 / (dr0 * (dr0 + dr1))
            drb = (dr1 - dr0) / (dr0 * dr1)
            drc = dr0 / (dr1 * (dr0 + dr1))

            y1_gradient_conj = dra * y1_conj_rm1 + drb * y1_conj + drc * y1_conj_rp1

        # Optimizations
        r2 = r * r
        y1y3_term_real = 2. * y1.real - llp1 * y3.real
        y1y3_term_imag = 2. * y1.imag - llp1 * y3.imag
        y1y3_term = y1y3_term_real + 1.0j * y1y3_term_imag
        y1y3_term_abs2 = (y1y3_term_real * y1y3_term_real) + (y1y3_term_imag * y1y3_term_imag)

        # Find bulk sensitivity
        if r == 0.:
            # The bulk sensitivity is not defined at r=0
            radial_sensitivity_to_bulk[r_i] = np.nan
        else:
            radial_sensitivity_to_bulk[r_i] = \
                (r2 / (np.abs(bulk + (4. / 3.) * shear)**2) *
                 np.abs(y2 - ((bulk - (2. / 3.) * shear) / r) * y1y3_term)**2) + \
                (2. * r * np.real(y1_gradient_conj * y1y3_term)) + \
                y1y3_term_abs2

    return radial_sensitivity_to_bulk