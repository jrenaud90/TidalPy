# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from libc.math cimport pi

from cpython.mem cimport PyMem_Malloc

# We need to use a custom csqrt function because Windows does not play nice with libc.complex cython library.
from TidalPy.utilities.math.complex cimport csqrt
from TidalPy.radial_solver.numerical.initial.common cimport takeuchi_phi_psi

cdef void takeuchi_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double bulk_modulus,
        double complex shear_modulus,
        unsigned char degree_l,
        double G_to_use,
        ssize_t num_ys, 
        double complex* initial_conditions
        ) noexcept nogil:
    """ Calculate the initial guess at the bottom of a solid layer using the dynamic assumption.

    This function uses the Takeuchi and Saito 1972 equations (Eq. 95-101).

    Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk and shear dissipation.

    References
    ----------
    TS72

    Parameters
    ----------
    radius : float
        Radius where the radial functions are calculated. [m]
    shear_modulus : Union[float, complex]
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : Union[float, complex]
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : float
        Density at `radius` [kg m-3]
    frequency : float
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    degree_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : np.ndarray
        The three independent solid guesses (sn1, sn2, sn3)

    """

    # Convert compressibility parameters
    cdef double complex lame
    lame = bulk_modulus - (2. / 3.) * shear_modulus

    # Constants (See Eqs. B13-B16 of KMN15)
    cdef double dynamic_term, gamma
    cdef double complex alpha2, beta2
    dynamic_term = frequency * frequency
    alpha2       = (lame + 2. * shear_modulus) / density
    beta2        = shear_modulus / density
    gamma        = 4. * pi * G_to_use * density / 3.

    # Optimizations
    cdef double r_inverse, r2, degree_l_dbl
    r_inverse    = 1. / radius
    r2           = radius * radius
    degree_l_dbl = <double> degree_l
    lp1          = degree_l_dbl + 1.
    lm1          = degree_l_dbl - 1.
    dlp1         = 2. * degree_l_dbl + 1.
    llp1         = degree_l_dbl * lp1
    dlp3         = 2. * degree_l_dbl + 3.

    # Helper functions
    cdef double complex k2_quad_pos, k2_quad_neg, k2_quad
    k2_quad_pos = (dynamic_term / beta2) + ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad_neg = (dynamic_term / beta2) - ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad = k2_quad_neg**2 + ((4. * llp1 * gamma**2) / (alpha2 * beta2))

    # TODO: TS74 has these flipped compared to KMN15. Going with TS74 for this func.
    cdef double complex k2_pos, k2_neg, k2_quad_sqrt
    k2_quad_sqrt = csqrt(k2_quad)
    k2_pos = (1. / 2.) * (k2_quad_pos - k2_quad_sqrt)
    k2_neg = (1. / 2.) * (k2_quad_pos + k2_quad_sqrt)

    cdef double complex f_k2_pos, f_k2_neg
    f_k2_pos = (beta2 * k2_pos - dynamic_term) / gamma
    f_k2_neg = (beta2 * k2_neg - dynamic_term) / gamma

    cdef double complex h_k2_pos, h_k2_neg
    h_k2_pos = f_k2_pos - lp1
    h_k2_neg = f_k2_neg - lp1

    # Calculate Takeuchi and Saito functions
    cdef double complex z_k2_pos, z_k2_neg
    z_k2_pos = k2_pos * r2
    z_k2_neg = k2_neg * r2

    cdef double complex phi_k2_pos, phi_lp1_k2_pos, psi_k2_pos
    cdef double complex phi_k2_neg, phi_lp1_k2_neg, psi_k2_neg
    takeuchi_phi_psi(z_k2_pos, degree_l, &phi_k2_pos, &phi_lp1_k2_pos, &psi_k2_pos)
    takeuchi_phi_psi(z_k2_neg, degree_l, &phi_k2_neg, &phi_lp1_k2_neg, &psi_k2_neg)

    # See Eq. 102 in TS72
    # y1, solutions 1--3
    initial_conditions[0 * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (.5 * degree_l_dbl * h_k2_pos * psi_k2_pos + f_k2_pos * phi_lp1_k2_pos)
    initial_conditions[1 * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (.5 * degree_l_dbl * h_k2_neg * psi_k2_neg + f_k2_neg * phi_lp1_k2_neg)
    initial_conditions[2 * num_ys + 0] = degree_l_dbl * radius**lm1
    # y2, solutions 1--3
    initial_conditions[0 * num_ys + 1] = \
        -(lame + 2. * shear_modulus) * radius**degree_l_dbl * f_k2_pos * phi_k2_pos + \
        (shear_modulus * radius**degree_l_dbl / dlp3) * (
            -degree_l_dbl * lm1 * h_k2_pos * psi_k2_pos + 2. * (2. * f_k2_pos + llp1) * phi_lp1_k2_pos
            )
    initial_conditions[1 * num_ys + 1] = \
        -(lame + 2. * shear_modulus) * radius**degree_l_dbl * f_k2_neg * phi_k2_neg + \
        (shear_modulus * radius**degree_l_dbl / dlp3) * (
            -degree_l_dbl * lm1 * h_k2_neg * psi_k2_neg + 2. * (2. * f_k2_neg + llp1) * phi_lp1_k2_neg
            )
    initial_conditions[2 * num_ys + 1] = \
        2. * shear_modulus * degree_l_dbl * lm1 * radius**(degree_l_dbl - 2.)
    # y3, solutions 1--3
    initial_conditions[0 * num_ys + 2] = \
        (-radius**lp1 / dlp3) * (0.5 * h_k2_pos * psi_k2_pos - phi_lp1_k2_pos)
    initial_conditions[1 * num_ys + 2] = \
        (-radius**lp1 / dlp3) * (0.5 * h_k2_neg * psi_k2_neg - phi_lp1_k2_neg)
    initial_conditions[2 * num_ys + 2] = \
        radius**lm1
    # y4, solutions 1--3
    initial_conditions[0 * num_ys + 3] = \
        shear_modulus * radius**degree_l_dbl * (
            phi_k2_pos - (1. / dlp3) * (lm1 * h_k2_pos * psi_k2_pos + 2. * (f_k2_pos + 1.) * phi_lp1_k2_pos)
            )
    initial_conditions[1 * num_ys + 3] = \
        shear_modulus * radius**degree_l_dbl * (
            phi_k2_neg - (1. / dlp3) * (lm1 * h_k2_neg * psi_k2_neg + 2. * (f_k2_neg + 1.) * phi_lp1_k2_neg)
            )
    initial_conditions[2 * num_ys + 3] = \
        2. * shear_modulus * lm1 * radius**(degree_l_dbl - 2.)
    # y5, solutions 1--3
    initial_conditions[0 * num_ys + 4] = \
        radius**(degree_l_dbl + 2.) * (
            (alpha2 * f_k2_pos - lp1 * beta2) / r2 - (3. * gamma * f_k2_pos / (2. * dlp3)) * psi_k2_pos
            )
    initial_conditions[1 * num_ys + 4] = \
        radius**(degree_l_dbl + 2.) * (
            (alpha2 * f_k2_neg - lp1 * beta2) / r2 - (3. * gamma * f_k2_neg / (2. * dlp3)) * psi_k2_neg
            )
    initial_conditions[2 * num_ys + 4] = \
        (degree_l_dbl * gamma - dynamic_term) * radius**degree_l_dbl
    # y6, solutions 1--3
    initial_conditions[0 * num_ys + 5] = \
        dlp1 * r_inverse * initial_conditions[0 * num_ys + 4] + \
        (3. * degree_l_dbl * gamma * h_k2_pos * radius**lp1 / (2. * dlp3)) * psi_k2_pos
    initial_conditions[1 * num_ys + 5] = \
        dlp1 * r_inverse * initial_conditions[1 * num_ys + 4] + \
        (3. * degree_l_dbl * gamma * h_k2_neg * radius**lp1 / (2. * dlp3)) * psi_k2_neg
    initial_conditions[2 * num_ys + 5] = \
        dlp1 * r_inverse * initial_conditions[2 * num_ys + 4] - \
        3. * degree_l_dbl * gamma * radius**lm1
