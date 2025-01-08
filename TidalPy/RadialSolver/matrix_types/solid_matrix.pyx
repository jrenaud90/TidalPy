# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
""" Fundamental Matrix and its inverse as defined in Sabadini, Vermeerson, & Cambiotti (2016) (Hereafter cited as SVC16)

Two versions are defined, one for a generic order-l and one for order-l=2 (for performance).

Assumptions
-----------
These matrices assume an incompressible body.

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
"""

from libcpp cimport bool as cpp_bool

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.constants cimport d_G, d_PI_DBL
from TidalPy.utilities.math.complex cimport cf_cinv


cdef void cf_fundamental_matrix(
        Py_ssize_t first_slice_index,
        Py_ssize_t num_radial_slices,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double complex* complex_shear_array_ptr,
        double complex* fundamental_mtx_ptr,
        double complex* inverse_fundamental_mtx_ptr,
        double complex* derivative_mtx_ptr,
        int degree_l = 2,
        double G_to_use = d_G
        ) noexcept nogil:
    """ Construct the fundamental matrix and its inverse for a generic order-l

    See Eq. 2.42 of SVC16

    Assumptions
    -----------
    - These matrices assume an incompressible body.

    Parameters
    ----------
    first_slice_index : Py_ssize_t
        Initial radial index to start populating matrices at.
    num_radial_slices : Py_ssize_t
        Number of radial slices
    radius_array_ptr : double
        Pointer to array of Radius values [m]
    density_array_ptr : double
        Pointer to array of Density at each radius [kg m-3]
    gravity_array_ptr : double
        Pointer to array_ptr of acceleration due to gravity at each radius [m s-2]
    complex_shear_array : double complex
        Pointer to array of Complex shear modulus at each radius [Pa]
    fundamental_mtx_ptr : double complex*
        _Return Value_
        6 x 6 matrix of double complex values
        Fundamental matrix used in the propagation technique
    inverse_fundamental_mtx_ptr : double complex*
        _Return Value_
        The inverse of the fundamental matrix used in the propagation technique
    derivative_mtx_ptr : double complex*
        _Return Value_
        The matrix, A, that satisfies the equation dy/dr = A * y
    degree_l : unsigned char, default=2
        Harmonic degree.
    G_to_use : double, default=d_G
        Gravitational constant used in calculations. Can be provided for non-dimensionalized solutions.
    """

    cdef Py_ssize_t slice_i, index_shift
    
    cdef double radius, gravity, density
    cdef double complex complex_shear
    
    cdef double r_inv, rl, rlp1, rlp2, rlp3, rnl, rnlm2, rlm1, rgp, piGp
    cdef double coeff, d_coeff_1, d_coeff_2, d_coeff_3, d_coeff_4, d_coeff_5, d_coeff_6
    cdef double complex mu_inv, rgp_s, r_s, pr_s

    # Degree-l Optimizations
    cdef double degree_l_dbl = <double> degree_l
    cdef double dlm1         = (2. * degree_l_dbl - 1.)
    cdef double l2p3lm1      = (degree_l_dbl**2 + 3. * degree_l_dbl - 1.)
    cdef double l2mlm3       = (degree_l_dbl**2 - degree_l_dbl - 3.)
    cdef double lp1          = (degree_l_dbl + 1.)
    cdef double lp2          = (degree_l_dbl + 2.)
    cdef double lp3          = (degree_l_dbl + 3.)
    cdef double l2m1         = (degree_l_dbl**2 - 1.)
    cdef double dlp1         = (2. * degree_l_dbl + 1.)
    cdef double dlp3         = (2. * degree_l_dbl + 3.)

    for slice_i in range(first_slice_index, num_radial_slices):

        # Shift index by 36 (for the inner 6x6 matrix)
        index_shift = slice_i * 36

        # Unpack radially dependent variables
        radius        = radius_array_ptr[slice_i]
        complex_shear = complex_shear_array_ptr[slice_i]
        gravity       = gravity_array_ptr[slice_i]
        density       = density_array_ptr[slice_i]

        # Radius-based optimizations
        r_inv  = 1. / radius
        mu_inv = cf_cinv(complex_shear)
        rl     = radius**degree_l
        rlp1   = radius**(degree_l + 1)
        rlp2   = radius**(degree_l + 2)
        rlp3   = radius**(degree_l + 3)
        rnl    = radius**(-degree_l)
        rnlm2  = radius**(-degree_l - 2)
        rlm1   = radius**(degree_l - 1)
        rgp    = radius * gravity * density
        rgp_s  = rgp * mu_inv
        r_s    = radius * mu_inv
        pr_s   = density * r_s
        piGp   = d_PI_DBL * G_to_use * density
        
        # D Coefficients
        coeff = (1. / dlp1)
        d_coeff_1 = coeff * lp1 / rlp1
        d_coeff_2 = coeff * degree_l_dbl * lp1 / (2. * dlm1 * rlm1)
        d_coeff_3 = coeff * 1. / rlm1
        d_coeff_4 = coeff * degree_l_dbl * rl
        d_coeff_5 = coeff * rlp2 * degree_l_dbl * lp1 / (2. * dlp3)
        d_coeff_6 = coeff * -rlp1

        
        # Build Fundamental Matrix (zeros do not need to be specifically stated as they were put in at initialization)
        #     Eq. 2.42 in SVC
        ## Row 1
        fundamental_mtx_ptr[index_shift + 0]  = degree_l_dbl * rlp1 / (2. * dlp3)
        fundamental_mtx_ptr[index_shift + 1]  = rlm1
        fundamental_mtx_ptr[index_shift + 2]  = 0.
        fundamental_mtx_ptr[index_shift + 3]  = lp1 * rnl / (2. * dlm1)
        fundamental_mtx_ptr[index_shift + 4]  = rnlm2
        fundamental_mtx_ptr[index_shift + 5]  = 0.

        ## Row 2
        fundamental_mtx_ptr[index_shift + 6]  = lp3 * rlp1 / (2. * dlp3 * lp1)
        fundamental_mtx_ptr[index_shift + 7]  = rlm1 / degree_l_dbl
        fundamental_mtx_ptr[index_shift + 8]  = 0.
        fundamental_mtx_ptr[index_shift + 9]  = (2. - degree_l_dbl) * rnl / (2. * degree_l_dbl * dlm1)
        fundamental_mtx_ptr[index_shift + 10] = -rnlm2 / lp1
        fundamental_mtx_ptr[index_shift + 11] = 0.

        ## Row 3
        # RECORD: Believe there is a typo in HH14, they have the radius^l only on one term instead of both.
        fundamental_mtx_ptr[index_shift + 12] = (degree_l_dbl * rgp + 2. * l2mlm3 * complex_shear) * rl / (2. * dlp3)
        fundamental_mtx_ptr[index_shift + 13] = (rgp + 2. * (degree_l_dbl - 1.) * complex_shear) * radius**(degree_l_dbl - 2.)
        fundamental_mtx_ptr[index_shift + 14] = -density * rl
        fundamental_mtx_ptr[index_shift + 15] = (lp1 * rgp - 2. * l2p3lm1 * complex_shear) / (2. * dlm1 * rlp1)
        fundamental_mtx_ptr[index_shift + 16] = (rgp - 2. * lp2 * complex_shear) / rlp3
        fundamental_mtx_ptr[index_shift + 17] = -density / rlp1

        ## Row 4
        fundamental_mtx_ptr[index_shift + 18] = degree_l_dbl * lp2 * complex_shear * rl / (dlp3 * lp1)
        fundamental_mtx_ptr[index_shift + 19] = 2 * (degree_l_dbl - 1.) * complex_shear * radius**(degree_l_dbl - 2.) / degree_l_dbl
        fundamental_mtx_ptr[index_shift + 20] = 0.
        fundamental_mtx_ptr[index_shift + 21] = l2m1 * complex_shear / (degree_l_dbl * dlm1 * rlp1)
        fundamental_mtx_ptr[index_shift + 22] = 2. * lp2 * complex_shear / (lp1 * rlp3)
        fundamental_mtx_ptr[index_shift + 23] = 0.

        ## Row 5
        fundamental_mtx_ptr[index_shift + 24] = 0.
        fundamental_mtx_ptr[index_shift + 25] = 0.
        fundamental_mtx_ptr[index_shift + 26] = -rl
        fundamental_mtx_ptr[index_shift + 27] = 0.
        fundamental_mtx_ptr[index_shift + 28] = 0.
        fundamental_mtx_ptr[index_shift + 29] = -1. / rlp1

        ## Row 6
        fundamental_mtx_ptr[index_shift + 30] = 2. * piGp * degree_l_dbl * rlp1 / dlp3
        fundamental_mtx_ptr[index_shift + 31] = 4. * piGp * rlm1
        fundamental_mtx_ptr[index_shift + 32] = -dlp1 * rlm1
        fundamental_mtx_ptr[index_shift + 33] = 2 * piGp * lp1 / (dlm1 * rl)
        fundamental_mtx_ptr[index_shift + 34] = 4. * piGp / rlp2
        fundamental_mtx_ptr[index_shift + 35] = 0.

        # Inverse of the Fundamental Matrix
        # This function manually defines the inverse matrix which is about 15--30% faster than a version that uses
        #    np.linalg.inv() to calculate the inverse matrix.
        #
        # From SVC16 Eq. 2.45: Fundamental Inverse = D_Mtx * Y^Bar_Mtx
        # D_Mtx is a diagonal matrix with
        # 1/(2l+1) * [ ... ]
        # We are going to multiple first and just write down the fundamental matrix inverse to avoid the additional D*Ybar
        #    calculation.

        ## Row 1
        inverse_fundamental_mtx_ptr[index_shift + 0]  = d_coeff_1 * (rgp_s - 2. * lp2)
        inverse_fundamental_mtx_ptr[index_shift + 1]  = d_coeff_1 * (2. * degree_l_dbl * lp2)
        inverse_fundamental_mtx_ptr[index_shift + 2]  = d_coeff_1 * (-r_s)
        inverse_fundamental_mtx_ptr[index_shift + 3]  = d_coeff_1 * (degree_l_dbl * r_s)
        inverse_fundamental_mtx_ptr[index_shift + 4]  = d_coeff_1 * (pr_s)
        inverse_fundamental_mtx_ptr[index_shift + 5]  = 0.

        ## Row 2
        inverse_fundamental_mtx_ptr[index_shift + 6]  = d_coeff_2 * (-rgp_s + 2. * l2p3lm1 / lp1)
        inverse_fundamental_mtx_ptr[index_shift + 7]  = d_coeff_2 * (-2. * l2m1)
        inverse_fundamental_mtx_ptr[index_shift + 8]  = d_coeff_2 * (r_s)
        inverse_fundamental_mtx_ptr[index_shift + 9]  = d_coeff_2 * ((2. - degree_l_dbl) * r_s)
        inverse_fundamental_mtx_ptr[index_shift + 10] = d_coeff_2 * (-pr_s)
        inverse_fundamental_mtx_ptr[index_shift + 11] = 0.

        ## Row 3
        inverse_fundamental_mtx_ptr[index_shift + 12] = d_coeff_3 * (4. * piGp)
        inverse_fundamental_mtx_ptr[index_shift + 13] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 14] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 15] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 16] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 17] = -d_coeff_3

        ## Row 4
        inverse_fundamental_mtx_ptr[index_shift + 18] = d_coeff_4 * (rgp_s + 2. * (degree_l_dbl - 1.))
        inverse_fundamental_mtx_ptr[index_shift + 19] = d_coeff_4 * (2. * l2m1)
        inverse_fundamental_mtx_ptr[index_shift + 20] = d_coeff_4 * (-r_s)
        inverse_fundamental_mtx_ptr[index_shift + 21] = d_coeff_4 * (-lp1 * r_s)
        inverse_fundamental_mtx_ptr[index_shift + 22] = d_coeff_4 * (pr_s)
        inverse_fundamental_mtx_ptr[index_shift + 23] = 0.

        ## Row 5
        inverse_fundamental_mtx_ptr[index_shift + 24] = d_coeff_5 * (-rgp_s - 2. * l2mlm3 / degree_l_dbl)
        inverse_fundamental_mtx_ptr[index_shift + 25] = d_coeff_5 * (-2. * degree_l_dbl * lp2)
        inverse_fundamental_mtx_ptr[index_shift + 26] = d_coeff_5 * (r_s)
        inverse_fundamental_mtx_ptr[index_shift + 27] = d_coeff_5 * (lp3 * r_s)
        inverse_fundamental_mtx_ptr[index_shift + 28] = d_coeff_5 * (-pr_s)
        inverse_fundamental_mtx_ptr[index_shift + 29] = 0.

        ## Row 6
        inverse_fundamental_mtx_ptr[index_shift + 30] = d_coeff_6 * (4. * piGp * radius)
        inverse_fundamental_mtx_ptr[index_shift + 31] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 32] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 33] = 0.
        inverse_fundamental_mtx_ptr[index_shift + 34] = d_coeff_6 * dlp1
        inverse_fundamental_mtx_ptr[index_shift + 35] = d_coeff_6 * (-radius)

        # Build derivative matrix
        # Defined in SV04 -- Only valid for the incompressible case.
        # See SVC16 Eq. 1.95
        #    Note: the lambda in SVC16 is defined as bulk_mod - (2. / 3.) * shear (Eq. 1.77; 2nd Lame parameter),
        #    for the incompressible assumption we will assume the ratio that SVC16 use (lambda / beta) -> 1 as K -> inf
        #    See SVC16 Eq. 1.95 for a compressible version. Take limit as K->inf to find below.
        ## Row 1
        derivative_mtx_ptr[index_shift + 0]  = -2. * r_inv
        derivative_mtx_ptr[index_shift + 1]  = degree_l_dbl * lp1 * r_inv
        derivative_mtx_ptr[index_shift + 2]  = 0.
        derivative_mtx_ptr[index_shift + 3]  = 0.
        derivative_mtx_ptr[index_shift + 4]  = 0.
        derivative_mtx_ptr[index_shift + 5]  = 0.

        ## Row 2
        derivative_mtx_ptr[index_shift + 6]  = -1. * r_inv
        derivative_mtx_ptr[index_shift + 7]  = r_inv
        derivative_mtx_ptr[index_shift + 8]  = 0.
        derivative_mtx_ptr[index_shift + 9]  = mu_inv
        derivative_mtx_ptr[index_shift + 10] = 0.
        derivative_mtx_ptr[index_shift + 11] = 0.

        ## Row 3
        derivative_mtx_ptr[index_shift + 12] = (4. * r_inv) * (3. * complex_shear * r_inv - density * gravity)
        derivative_mtx_ptr[index_shift + 13] = -degree_l_dbl * lp1 * r_inv * (6. * complex_shear * r_inv - density * gravity)
        derivative_mtx_ptr[index_shift + 14] = 0.
        derivative_mtx_ptr[index_shift + 15] = degree_l_dbl * lp1 * r_inv
        derivative_mtx_ptr[index_shift + 16] = -density * lp1 * r_inv
        derivative_mtx_ptr[index_shift + 17] = density

        ## Row 4
        derivative_mtx_ptr[index_shift + 18] = (-1. * r_inv) * (6. * complex_shear * r_inv - density * gravity)
        derivative_mtx_ptr[index_shift + 19] = 2. * (2. * degree_l_dbl**2 + 2. * degree_l_dbl - 1.) * complex_shear * (r_inv**2)
        derivative_mtx_ptr[index_shift + 20] = -r_inv
        derivative_mtx_ptr[index_shift + 21] = -3. * r_inv
        derivative_mtx_ptr[index_shift + 22] = density * r_inv
        derivative_mtx_ptr[index_shift + 23] = 0.

        ## Row 5
        derivative_mtx_ptr[index_shift + 24] = -4. * piGp
        derivative_mtx_ptr[index_shift + 25] = 0.
        derivative_mtx_ptr[index_shift + 26] = 0.
        derivative_mtx_ptr[index_shift + 27] = 0.
        derivative_mtx_ptr[index_shift + 28] = -lp1 * r_inv
        derivative_mtx_ptr[index_shift + 29] = 1.

        ## Row 6
        derivative_mtx_ptr[index_shift + 30] = -4. * piGp * (degree_l_dbl + 1) * r_inv
        derivative_mtx_ptr[index_shift + 31] = 4. * piGp * degree_l_dbl * lp1 * r_inv
        derivative_mtx_ptr[index_shift + 32] = 0.
        derivative_mtx_ptr[index_shift + 33] = 0.
        derivative_mtx_ptr[index_shift + 34] = 0.
        derivative_mtx_ptr[index_shift + 35] = (degree_l_dbl - 1.) * r_inv


def fundamental_matrix(
        double[::1] radius_array_view,
        double[::1] density_array_view,
        double[::1] gravity_array_view,
        double complex[::1] complex_shear_array_view,
        int degree_l = 2,
        double G_to_use = d_G,
        cpp_bool perform_checks = True
        ):
    """ Construct the fundamental matrix and its inverse using harmonic degree l.

    See Eq. 2.42 of SVC16

    Assumptions
    -----------
    - These matrices assume an incompressible body.

    Parameters
    ----------
    radius_array_view : double
        Pointer to array of Radius values [m]
    density_array_view : double
        Pointer to array of Density at each radius [kg m-3]
    gravity_array_view : double
        Pointer to array of acceleration due to gravity at each radius [m s-2]
    complex_shear_array_view : double complex
        Pointer to array of Complex shear modulus at each radius [Pa]
    degree_l : int, default=2
        Harmonic degree.
    G_to_use : double, default=d_G
        Gravitational constant used in calculations. Can be provided for non-dimensionalized solutions.
    perform_checks : bool, default=True
        If True, then checks will be performed on the input arguments to check for issues.
    """ 

    cdef Py_ssize_t num_radial_slices
    num_radial_slices = len(radius_array_view)

    # Check for unexpected shapes and sizes of return matrices
    if perform_checks:
        if len(density_array_view) != num_radial_slices:
            raise ValueError('Unexpected size encountered for density array.')
        if len(gravity_array_view) != num_radial_slices:
            raise ValueError('Unexpected size encountered for gravity array.')
        if len(complex_shear_array_view) != num_radial_slices:
            raise ValueError('Unexpected size encountered for complex shear array.')
    
    # Build output arrays
    cdef cnp.ndarray[cnp.complex128_t, ndim=3] fundamental_mtx_arr         = np.empty((num_radial_slices, 6, 6), dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=3] inverse_fundamental_mtx_arr = np.empty((num_radial_slices, 6, 6), dtype=np.complex128, order='C')
    cdef cnp.ndarray[cnp.complex128_t, ndim=3] derivative_mtx_arr          = np.empty((num_radial_slices, 6, 6), dtype=np.complex128, order='C')

    cdef double complex[:, :, ::1] fundamental_mtx_view         = fundamental_mtx_arr
    cdef double complex[:, :, ::1] inverse_fundamental_mtx_view = inverse_fundamental_mtx_arr
    cdef double complex[:, :, ::1] derivative_mtx_view          = derivative_mtx_arr

    # Call cdef function.
    cf_fundamental_matrix(
        0,
        num_radial_slices,
        &radius_array_view[0],
        &density_array_view[0],
        &gravity_array_view[0],
        &complex_shear_array_view[0],
        &fundamental_mtx_view[0, 0, 0],
        &inverse_fundamental_mtx_view[0, 0, 0],
        &derivative_mtx_view[0, 0, 0],
        degree_l,
        G_to_use
    )

    return fundamental_mtx_arr, inverse_fundamental_mtx_arr, derivative_mtx_arr
