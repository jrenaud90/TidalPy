# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from libc.math cimport pi

# We need to use a custom cf_csqrt function because Windows does not play nice with libc.complex cython library.
from TidalPy.utilities.math.complex cimport cf_csqrt
from TidalPy.RadialSolver.starting.common cimport cf_takeuchi_phi_psi

########################################################################################################################
#### Solid Layers
########################################################################################################################


cdef void cf_takeuchi_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double bulk_modulus,
        double complex shear_modulus,
        unsigned char degree_l,
        double G_to_use,
        ssize_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """ Calculate the starting guess at the bottom of a solid layer using the dynamic assumption.

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
    cdef double r_inverse, r2, degree_l_dbl, lp1, lm1, dlp1, llp1, dlp3
    r_inverse    = 1. / radius
    r2           = radius * radius
    degree_l_dbl = <double> degree_l
    lp1          = degree_l_dbl + 1.
    lm1          = degree_l_dbl - 1.
    dlp1         = 2. * degree_l_dbl + 1.
    dlp3         = 2. * degree_l_dbl + 3.
    llp1         = degree_l_dbl * lp1

    # Helper functions
    cdef double complex k2_quad_pos, k2_quad_neg, k2_quad, k2_quad_sqrt
    k2_quad_pos  = (dynamic_term / beta2) + ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad_neg  = (dynamic_term / beta2) - ((dynamic_term + 4. * gamma) / alpha2)
    k2_quad      = (k2_quad_neg * k2_quad_neg) + ((4. * llp1 * gamma**2) / (alpha2 * beta2))
    k2_quad_sqrt = cf_csqrt(k2_quad)

    # QUESTION: (Issue #43) TS74 has these flipped compared to KMN15. Going with TS74 for this func.
    # See the -/+ order in TS72 EQ. 99
    cdef size_t neg_index = 0
    cdef size_t pos_index = 1
    cdef double complex k2_pos, k2_neg
    k2_neg = (1. / 2.) * (k2_quad_pos - k2_quad_sqrt)
    k2_pos = (1. / 2.) * (k2_quad_pos + k2_quad_sqrt)

    cdef double complex f_pos, f_neg
    f_pos = (beta2 * k2_pos - dynamic_term) / gamma
    f_neg = (beta2 * k2_neg - dynamic_term) / gamma

    cdef double complex h_pos, h_neg
    h_pos = f_pos - lp1
    h_neg = f_neg - lp1

    # Calculate Takeuchi and Saito functions
    cdef double complex z2_pos, z2_neg
    z2_pos = k2_pos * r2
    z2_neg = k2_neg * r2

    cdef double complex phi_pos, phi_lp1_pos, psi_pos
    cdef double complex phi_neg, phi_lp1_neg, psi_neg
    cf_takeuchi_phi_psi(z2_pos, degree_l, &phi_pos, &phi_lp1_pos, &psi_pos)
    cf_takeuchi_phi_psi(z2_neg, degree_l, &phi_neg, &phi_lp1_neg, &psi_neg)

    # See Eq. 102 in TS72
    
    # y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (.5 * degree_l_dbl * h_pos * psi_pos + f_pos * phi_lp1_pos)
    starting_conditions_ptr[neg_index * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (.5 * degree_l_dbl * h_neg * psi_neg + f_neg * phi_lp1_neg)
    starting_conditions_ptr[2 * num_ys + 0] = degree_l_dbl * radius**lm1
    
    # y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] = \
        -(lame + 2. * shear_modulus) * radius**degree_l_dbl * f_pos * phi_pos + \
        (shear_modulus * radius**degree_l_dbl / dlp3) * (
            -degree_l_dbl * lm1 * h_pos * psi_pos + 2. * (2. * f_pos + llp1) * phi_lp1_pos
            )
    starting_conditions_ptr[neg_index * num_ys + 1] = \
        -(lame + 2. * shear_modulus) * radius**degree_l_dbl * f_neg * phi_neg + \
        (shear_modulus * radius**degree_l_dbl / dlp3) * (
            -degree_l_dbl * lm1 * h_neg * psi_neg + 2. * (2. * f_neg + llp1) * phi_lp1_neg
            )
    starting_conditions_ptr[2 * num_ys + 1] = \
        2. * shear_modulus * degree_l_dbl * lm1 * radius**(degree_l_dbl - 2.)
    
    # y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] = \
        (-radius**lp1 / dlp3) * (0.5 * h_pos * psi_pos - phi_lp1_pos)
    starting_conditions_ptr[neg_index * num_ys + 2] = \
        (-radius**lp1 / dlp3) * (0.5 * h_neg * psi_neg - phi_lp1_neg)
    starting_conditions_ptr[2 * num_ys + 2] = \
        radius**lm1
    
    # y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] = \
        shear_modulus * radius**degree_l_dbl * (
            phi_pos - (1. / dlp3) * (lm1 * h_pos * psi_pos + 2. * (f_pos + 1.) * phi_lp1_pos)
            )
    starting_conditions_ptr[neg_index * num_ys + 3] = \
        shear_modulus * radius**degree_l_dbl * (
            phi_neg - (1. / dlp3) * (lm1 * h_neg * psi_neg + 2. * (f_neg + 1.) * phi_lp1_neg)
            )
    starting_conditions_ptr[2 * num_ys + 3] = \
        2. * shear_modulus * lm1 * radius**(degree_l_dbl - 2.)
    
    # y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] = \
        radius**(degree_l_dbl + 2.) * (
            (alpha2 * f_pos - lp1 * beta2) / r2 - (3. * gamma * f_pos / (2. * dlp3)) * psi_pos
            )
    starting_conditions_ptr[neg_index * num_ys + 4] = \
        radius**(degree_l_dbl + 2.) * (
            (alpha2 * f_neg - lp1 * beta2) / r2 - (3. * gamma * f_neg / (2. * dlp3)) * psi_neg
            )
    starting_conditions_ptr[2 * num_ys + 4] = \
        (degree_l_dbl * gamma - dynamic_term) * radius**degree_l_dbl
    
    # y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] = \
        dlp1 * r_inverse * starting_conditions_ptr[0 * num_ys + 4] + \
        (3. * degree_l_dbl * gamma * h_pos * radius**lp1 / (2. * dlp3)) * psi_pos
    starting_conditions_ptr[neg_index * num_ys + 5] = \
        dlp1 * r_inverse * starting_conditions_ptr[1 * num_ys + 4] + \
        (3. * degree_l_dbl * gamma * h_neg * radius**lp1 / (2. * dlp3)) * psi_neg
    starting_conditions_ptr[2 * num_ys + 5] = \
        dlp1 * r_inverse * starting_conditions_ptr[2 * num_ys + 4] - \
        3. * degree_l_dbl * gamma * radius**lm1


cdef void cf_takeuchi_solid_static_compressible(
        double radius,
        double density,
        double bulk_modulus,
        double complex shear_modulus,
        unsigned char degree_l,
        double G_to_use,
        ssize_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """ Calculate the starting guess at the bottom of a solid layer using the static assumption.

    This function uses the Takeuchi and Saito 1972 equations (Eq. 95-101).

    Using the static assumption in a solid layer results in two independent solutions for the radial derivative.

    These independent solutions allow for a general tidal harmonic l, for static tides (w = 0), compressibility, and
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
    cdef double complex alpha2, beta2
    cdef double gamma
    alpha2 = (lame + 2. * shear_modulus) / density
    beta2  = shear_modulus / density
    gamma   = 4. * pi * G_to_use * density / 3.

    # Optimizations
    cdef double r_inverse, r2, degree_l_dbl, lp1, lm1, dlp1, llp1, dlp3
    r_inverse    = 1. / radius
    r2           = radius * radius
    degree_l_dbl = <double> degree_l
    lp1          = degree_l_dbl + 1.
    lm1          = degree_l_dbl - 1.
    dlp1         = 2. * degree_l_dbl + 1.
    dlp3         = 2. * degree_l_dbl + 3.
    llp1         = degree_l_dbl * lp1

    # Helper functions
    cdef double complex k2_quad_pos, k2_quad, k2_quad_sqrt
    k2_quad_pos  = (4. * gamma / alpha2)
    k2_quad      = k2_quad_pos**2 + ((4. * llp1 * gamma**2) / (alpha2 * beta2))
    k2_quad_sqrt = cf_csqrt(k2_quad)

    # QUESTION: (Issue #43) TS74 has these flipped compared to KMN15. Going with TS74 for this func.
    # See the -/+ order in TS72 EQ. 99
    cdef size_t neg_index = 0
    cdef size_t pos_index = 1
    cdef double complex k2_pos, k2_neg
    k2_neg = (1. / 2.) * (k2_quad_pos - k2_quad_sqrt)
    k2_pos = (1. / 2.) * (k2_quad_pos + k2_quad_sqrt)

    cdef double complex f_pos, f_neg
    f_pos = (beta2 * k2_pos) / gamma
    f_neg = (beta2 * k2_neg) / gamma

    cdef double complex h_pos, h_neg
    h_pos = f_pos - lp1
    h_neg = f_neg - lp1

    # Calculate Takeuchi and Saito functions
    cdef double complex z2_pos, z2_neg
    z2_pos = k2_pos * r2
    z2_neg = k2_neg * r2

    cdef double complex phi_pos, phi_lp1_pos, psi_pos
    cdef double complex phi_neg, phi_lp1_neg, psi_neg
    cf_takeuchi_phi_psi(z2_pos, degree_l, &phi_pos, &phi_lp1_pos, &psi_pos)
    cf_takeuchi_phi_psi(z2_neg, degree_l, &phi_neg, &phi_lp1_neg, &psi_neg)

    # See Eq. 102 in TS72
    # y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (0.5 * degree_l * h_pos * psi_pos + f_pos * phi_lp1_pos)
    starting_conditions_ptr[neg_index * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (0.5 * degree_l * h_neg * psi_neg + f_neg * phi_lp1_neg)
    starting_conditions_ptr[2 * num_ys + 0] = \
        degree_l * radius**lm1

    # y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] = \
        -(lame + 2. * shear_modulus) * radius**degree_l * f_pos * phi_pos + \
        (shear_modulus * radius**degree_l / dlp3) * (
            -degree_l * lm1 * h_pos * psi_pos + 2. * (2. * f_pos + llp1) * phi_lp1_pos
            )
    starting_conditions_ptr[neg_index * num_ys + 1] = \
        -(lame + 2. * shear_modulus) * radius**degree_l * f_neg * phi_neg + \
        (shear_modulus * radius**degree_l / dlp3) * (
            -degree_l * lm1 * h_neg * psi_neg + 2. * (2. * f_neg + llp1) * phi_lp1_neg
            )
    starting_conditions_ptr[2 * num_ys + 1] = \
        2. * shear_modulus * degree_l * lm1 * radius**(degree_l - 2.)

    # y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] = \
        (-radius**lp1 / dlp3) * (0.5 * h_pos * psi_pos - phi_lp1_pos)
    starting_conditions_ptr[neg_index * num_ys + 2] = \
        (-radius**lp1 / dlp3) * (0.5 * h_neg * psi_neg - phi_lp1_neg)
    starting_conditions_ptr[2 * num_ys + 2] = \
        radius**lm1

    # y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] = \
        shear_modulus * radius**degree_l * \
        (phi_pos - (1. / dlp3) * (lm1 * h_pos * psi_pos + 2. * (f_pos + 1.) * phi_lp1_pos))
    starting_conditions_ptr[neg_index * num_ys + 3] = \
        shear_modulus * radius**degree_l * \
        (phi_neg - (1. / dlp3) * (lm1 * h_neg * psi_neg + 2. * (f_neg + 1.) * phi_lp1_neg))
    starting_conditions_ptr[2 * num_ys + 3] = \
        2. * shear_modulus * lm1 * radius**(degree_l - 2.)

    # y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] = \
        radius**(degree_l + 2.) * \
        ((alpha2 * f_pos - lp1 * beta2) / r2 - (3. * gamma * f_pos / (2. * dlp3)) * psi_pos)
    starting_conditions_ptr[neg_index * num_ys + 4] = \
        radius**(degree_l + 2.) * \
        ((alpha2 * f_neg - lp1 * beta2) / r2 - (3. * gamma * f_neg / (2. * dlp3)) * psi_neg)
    starting_conditions_ptr[2 * num_ys + 4] = \
        degree_l * gamma * radius**degree_l

    # y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] = \
        dlp1 * r_inverse * starting_conditions_ptr[0 * num_ys + 4] + \
        (3. * degree_l * gamma * h_pos * radius**lp1 / (2. * dlp3)) * psi_pos
    starting_conditions_ptr[neg_index * num_ys + 5] = \
        dlp1 * r_inverse * starting_conditions_ptr[1 * num_ys + 4] + \
        (3. * degree_l * gamma * h_neg * radius**lp1 / (2. * dlp3)) * psi_neg
    starting_conditions_ptr[2 * num_ys + 5] = \
        dlp1 * r_inverse * starting_conditions_ptr[2 * num_ys + 4] - \
        3. * degree_l * gamma * radius**lm1


########################################################################################################################
#### Liquid Layers
########################################################################################################################


cdef void cf_takeuchi_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double bulk_modulus,
        unsigned char degree_l,
        double G_to_use,
        ssize_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """ Calculate the starting guess at the bottom of a liquid layer using the dynamic assumption.

    This function uses the Takeuchi and Saito 1972 equations (Eq. 95-101).

    Using the dynamic assumption in a liquid layer results in two independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk dissipation (no shear dissipation within liquid layers).

    References
    ----------
    TS72

    Parameters
    ----------
    radius : float
        Radius where the radial functions are calculated. [m]
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
        The two independent liquid guesses (sn1, sn2)

    """

    # Convert compressibility parameters
    cdef double lame
    lame = bulk_modulus

    # Constants (See Eqs. B13-B16 of KMN15)
    cdef double dynamic_term, alpha2, gamma
    dynamic_term = frequency * frequency
    alpha2       = lame / density
    gamma        = 4. * pi * G_to_use * density / 3.

    # Optimizations
    cdef double r_inverse, r2, degree_l_dbl, lp1, lm1, dlp1, llp1, dlp3
    r_inverse    = 1. / radius
    r2           = radius * radius
    degree_l_dbl = <double> degree_l
    lp1          = degree_l_dbl + 1.
    lm1          = degree_l_dbl - 1.
    dlp1         = 2. * degree_l_dbl + 1.
    dlp3         = 2. * degree_l_dbl + 3.
    llp1         = degree_l_dbl * lp1

    # Helper functions
    cdef double k2, f, h
    k2 = (1. / alpha2) * (dynamic_term + 4. * gamma - llp1 * gamma**2 / dynamic_term)
    f  = -dynamic_term / gamma
    h  = f - lp1

    # Calculate Takeuchi and Saito functions
    cdef double complex z, phi, phi_lp1, psi
    z = k2 * r2
    cf_takeuchi_phi_psi(z, degree_l, &phi, &phi_lp1, &psi)

    # Found by setting mu=0 in Eq. 102 of TS72
    # y1, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 0] = \
        ((-radius**lp1) / dlp3) * (0.5 * degree_l * h * psi + f * phi_lp1)
    starting_conditions_ptr[1 * num_ys + 0] = \
        degree_l * radius**lm1

    # y2, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 1] = \
        -lame * radius**degree_l * f * phi
    starting_conditions_ptr[1 * num_ys + 1] = \
        0.0

    # y5, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 2] = \
        radius**(degree_l + 2.) * ((alpha2 * f / r2) - (3. * gamma * f / (2. * dlp3)) * psi)
    starting_conditions_ptr[1 * num_ys + 2] = \
        (degree_l * gamma - dynamic_term) * radius**degree_l

    # y6, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 3] = \
        dlp1 * r_inverse * starting_conditions_ptr[0 * num_ys + 2] + (3. * degree_l * gamma * h * radius**lp1 / (2. * dlp3)) * psi
    starting_conditions_ptr[1 * num_ys + 3] = \
        dlp1 * r_inverse * starting_conditions_ptr[1 * num_ys + 2] - 3. * degree_l * gamma * radius**lm1
