# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
""" starting conditions at the center of the planet based off of Kamata et al. (2015). """

from TidalPy.constants cimport d_PI_DBL
from TidalPy.utilities.math.complex cimport cf_csqrt
from TidalPy.RadialSolver.starting.common cimport cf_z_calc

########################################################################################################################
#### Solid Layers
########################################################################################################################


cdef void cf_kamata_solid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """ Calculate the starting guess at the bottom of a solid layer using the dynamic assumption.

    This function uses the Kamata et al (2015; JGR:P) equations (Eq. B1-B16).

    Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk and shear dissipation.

    References
    ----------
    KMN15 Eqs. B1-B16

    Parameters
    ----------
    frequency : double
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    radius : double
        Radius where the radial functions are calculated. [m]
    density : double
        Density at `radius` [kg m-3]
    bulk_modulus : double complex
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    shear_modulus : double complex
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    degree_l : unsigned char
        Tidal harmonic order.
    G_to_use : double
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    num_ys : ssize_t
        Number of radial solutions for this layer type.
    starting_conditions_ptr : double complex*, <Output>
        Desired starting conditions for this layer.
        Three independent solid guesses (sn1, sn2, sn3)
    """

    # Convert compressibility parameters
    cdef double complex lame = bulk_modulus - (2. / 3.) * shear_modulus

    # Constants
    cdef double gamma          = 4. * d_PI_DBL * G_to_use * density / 3.
    cdef double dynamic_term   = frequency * frequency
    cdef double complex alpha2 = (lame + 2. * shear_modulus) / density
    cdef double complex beta2  = shear_modulus / density

    # Optimizations
    cdef double r_inverse    = 1. / radius
    cdef double r2_inverse   = r_inverse * r_inverse
    cdef double r2           = radius * radius
    cdef double degree_l_dbl = <double> degree_l
    cdef double lp1          = degree_l_dbl + 1.0
    cdef double dlp1         = 2.0 * degree_l_dbl + 1.0
    cdef double llp1         = degree_l_dbl * lp1

    # Helper functions
    cdef double complex k2_quad_pos = (dynamic_term / beta2) + ((dynamic_term + 4. * gamma) / alpha2)
    cdef double complex k2_quad_neg = (dynamic_term / beta2) - ((dynamic_term + 4. * gamma) / alpha2)
    cdef double complex k2_quad     = (k2_quad_neg * k2_quad_neg) + \
        ((4. * degree_l * (degree_l + 1) * (gamma * gamma)) / (alpha2 * beta2))

    # QUESTION: (Issue #43) KMN15 has these flipped compared to TS72. Going with  KMN15 for this func.
    cdef size_t neg_index = 1
    cdef size_t pos_index = 0
    cdef double complex k2_quad_sqrt = cf_csqrt(k2_quad)
    cdef double complex k2_pos = (1. / 2.) * (k2_quad_pos + k2_quad_sqrt)
    cdef double complex k2_neg = (1. / 2.) * (k2_quad_pos - k2_quad_sqrt)

    cdef double complex f_k2_pos = (beta2 * k2_pos - dynamic_term) / gamma
    cdef double complex f_k2_neg = (beta2 * k2_neg - dynamic_term) / gamma

    cdef double complex h_k2_pos = f_k2_pos - lp1
    cdef double complex h_k2_neg = f_k2_neg - lp1

    cdef double complex z_k2_pos = cf_z_calc(k2_pos * r2, degree_l=degree_l)
    cdef double complex z_k2_neg = cf_z_calc(k2_neg * r2, degree_l=degree_l)

    # See Eqs. B1-B12 of KMN15
    
    # y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] = \
        -f_k2_pos * z_k2_pos * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 0] = \
        -f_k2_neg * z_k2_neg * r_inverse
    starting_conditions_ptr[2 * num_ys + 0] = \
        degree_l_dbl * r_inverse
    
    # y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] = \
        -density * f_k2_pos * alpha2 * k2_pos + (2. * shear_modulus * r2_inverse) * (2. * f_k2_pos + llp1) * z_k2_pos
    starting_conditions_ptr[neg_index * num_ys + 1] = \
        -density * f_k2_neg * alpha2 * k2_neg + (2. * shear_modulus * r2_inverse) * (2. * f_k2_neg + llp1) * z_k2_neg
    starting_conditions_ptr[2 * num_ys + 1] = \
        2. * shear_modulus * degree_l_dbl * (degree_l_dbl - 1) * r2_inverse
    
    # y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] = \
        z_k2_pos * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 2] = \
        z_k2_neg * r_inverse
    starting_conditions_ptr[2 * num_ys + 2] = \
        r_inverse
    
    # y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] = \
        shear_modulus * k2_pos - (2. * shear_modulus * r2_inverse) * (f_k2_pos + 1.) * z_k2_pos
    starting_conditions_ptr[neg_index * num_ys + 3] = \
        shear_modulus * k2_neg - (2. * shear_modulus * r2_inverse) * (f_k2_neg + 1.) * z_k2_neg
    starting_conditions_ptr[2 * num_ys + 3] = \
        2. * shear_modulus * (degree_l_dbl - 1.) * r2_inverse
    
    # y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] = \
        3. * gamma * f_k2_pos - h_k2_pos * (degree_l_dbl * gamma - dynamic_term)
    starting_conditions_ptr[neg_index * num_ys + 4] = \
        3. * gamma * f_k2_neg - h_k2_neg * (degree_l_dbl * gamma - dynamic_term)
    starting_conditions_ptr[2 * num_ys + 4] = \
        degree_l * gamma - dynamic_term
    
    # y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[pos_index * num_ys + 4] * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[neg_index * num_ys + 4] * r_inverse
    starting_conditions_ptr[2 * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[2 * num_ys + 4] * r_inverse - (3. * degree_l_dbl * gamma * r_inverse)


cdef void cf_kamata_solid_static_compressible(
        double radius,
        double density,
        double complex bulk_modulus,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """ Calculate the starting guess at the bottom of a solid layer using the static assumption.

    This function uses the Kamata et al. (2015; JGR:P) equations (Eq. B1-B16).

    Using the static assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for static tides (w = 0), compressibility, and
       bulk and shear dissipation.

    References
    ----------
    KMN15

    Parameters
    ----------
    radius : double
        Radius where the radial functions are calculated. [m]
    density : double
        Density at `radius` [kg m-3]
    bulk_modulus : double complex
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    shear_modulus : double complex
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    degree_l : unsigned char
        Tidal harmonic order.
    G_to_use : double
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    num_ys : ssize_t
        Number of radial solutions for this layer type.
    starting_conditions_ptr : double complex*, <Output>
        Desired starting conditions for this layer.
        Three independent solid guesses (sn1, sn2, sn3)
    """

    # Convert compressibility parameters
    cdef double complex lame = bulk_modulus - (2. / 3.) * shear_modulus

    # Constants
    cdef double gamma          = 4. * d_PI_DBL * G_to_use * density / 3.
    cdef double complex alpha2 = (lame + 2. * shear_modulus) / density
    cdef double complex beta2  = shear_modulus / density

    # Optimizations
    cdef double r_inverse    = 1. / radius
    cdef double r2_inverse   = r_inverse * r_inverse
    cdef double r2           = radius * radius
    cdef double degree_l_dbl = <double>degree_l
    cdef double lp1          = degree_l_dbl + 1.
    cdef double dlp1         = 2.0 * degree_l_dbl + 1.0
    cdef double llp1         = degree_l_dbl * lp1

    # Helper functions
    cdef double complex k2_quad_pos = 4. * gamma / alpha2
    cdef double complex k2_quad_neg = -k2_quad_pos
    cdef double complex k2_quad = k2_quad_neg**2 + ((4. * degree_l * lp1 * gamma**2) / (alpha2 * beta2))

    # QUESTION: (Issue #43) KMN15 has these flipped compared to TS72. Going with  KMN15 for this func.
    cdef size_t neg_index = 1
    cdef size_t pos_index = 0
    cdef double complex k2_quad_sqrt = cf_csqrt(k2_quad)
    cdef double complex k2_pos = (1. / 2.) * (k2_quad_pos + k2_quad_sqrt)
    cdef double complex k2_neg = (1. / 2.) * (k2_quad_pos - k2_quad_sqrt)
    
    cdef double complex f_k2_pos = beta2 * k2_pos / gamma
    cdef double complex f_k2_neg = beta2 * k2_neg / gamma

    cdef double complex h_k2_pos = f_k2_pos - lp1
    cdef double complex h_k2_neg = f_k2_neg - lp1

    cdef double complex z_k2_pos = cf_z_calc(k2_pos * r2, degree_l=degree_l)
    cdef double complex z_k2_neg = cf_z_calc(k2_neg * r2, degree_l=degree_l)

    # See Eqs. B1-B12 of KMN15
    
    # y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] = \
        -f_k2_pos * z_k2_pos * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 0] = \
        -f_k2_neg * z_k2_neg * r_inverse
    starting_conditions_ptr[2 * num_ys + 0] = \
        degree_l * r_inverse
    
    # y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] = \
        -density * f_k2_pos * alpha2 * k2_pos + (2. * shear_modulus * r2_inverse) * (2. * f_k2_pos + llp1) * z_k2_pos
    starting_conditions_ptr[neg_index * num_ys + 1] = \
        -density * f_k2_neg * alpha2 * k2_neg + (2. * shear_modulus * r2_inverse) * (2. * f_k2_neg + llp1) * z_k2_neg
    starting_conditions_ptr[2 * num_ys + 1] = \
        2. * shear_modulus * degree_l * (degree_l - 1) * r2_inverse
    
    # y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] = \
        z_k2_pos * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 2] = \
        z_k2_neg * r_inverse
    starting_conditions_ptr[2 * num_ys + 2] = \
        r_inverse
    
    # y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] = \
        shear_modulus * k2_pos - (2. * shear_modulus * r2_inverse) * (f_k2_pos + 1.) * z_k2_pos
    starting_conditions_ptr[neg_index * num_ys + 3] = \
        shear_modulus * k2_neg - (2. * shear_modulus * r2_inverse) * (f_k2_neg + 1.) * z_k2_neg
    starting_conditions_ptr[2 * num_ys + 3] = \
        2. * shear_modulus * (degree_l - 1.) * r2_inverse
    
    # y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] = \
        3. * gamma * f_k2_pos - h_k2_pos * (degree_l * gamma)
    starting_conditions_ptr[neg_index * num_ys + 4] = \
        3. * gamma * f_k2_neg - h_k2_neg * (degree_l * gamma)
    starting_conditions_ptr[2 * num_ys + 4] = \
        degree_l * gamma
    
    # y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[pos_index * num_ys + 4] * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[neg_index * num_ys + 4] * r_inverse
    starting_conditions_ptr[2 * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[2 * num_ys + 4] * r_inverse - (3. * degree_l * gamma * r_inverse)


cdef void cf_kamata_solid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        double complex shear_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """ Calculate the starting guess at the bottom of a solid layer using the dynamic and incompressible assumption.

    This function uses the Kamata et al. (2015; JGR:P) equations (Eq. B17-B28).

    Using the dynamic assumption in a solid layer results in three independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), incompressible, and
       bulk and shear dissipation.

    References
    ----------
    KMN15 Eqs. B17-28

    Parameters
    ----------
    frequency : double
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    radius : double
        Radius where the radial functions are calculated. [m]
    density : double
        Density at `radius` [kg m-3]
    shear_modulus : double complex
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    degree_l : unsigned char
        Tidal harmonic order.
    G_to_use : double
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    num_ys : ssize_t
        Number of radial solutions for this layer type.
    starting_conditions_ptr : double complex*, <Output>
        Desired starting conditions for this layer.
        Three independent solid guesses (sn1, sn2, sn3)

    """

    # Constants
    cdef double gamma         = 4. * d_PI_DBL * G_to_use * density / 3.
    cdef double dynamic_term  = frequency * frequency
    cdef double complex beta2 = shear_modulus / density

    # Optimizations
    cdef double r_inverse    = 1. / radius
    cdef double r2_inverse   = r_inverse * r_inverse
    cdef double r2           = radius * radius
    cdef double degree_l_dbl = <double>degree_l
    cdef double lp1          = degree_l_dbl + 1.
    cdef double lm1          = degree_l_dbl - 1.
    cdef double dlp1         = 2.0 * degree_l_dbl + 1.0
    cdef double llp1         = degree_l_dbl * lp1

    # QUESTION: (Issue #43) KMN15 has these flipped compared to TS72. Going with  KMN15 for this func.
    cdef size_t neg_index = 1
    cdef size_t pos_index = 0
    cdef double complex k2_pos   = dynamic_term / beta2
    cdef double complex f_k2_neg = -dynamic_term / gamma
    cdef double complex h_k2_neg = f_k2_neg - lp1
    cdef double complex z_k2_pos = cf_z_calc(k2_pos * r2, degree_l=degree_l)

    # See Eqs. B17-B28 of KMN15
    
    # y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] = \
        0.
    starting_conditions_ptr[neg_index * num_ys + 0] = \
        0.
    starting_conditions_ptr[2 * num_ys + 0] = \
        degree_l_dbl * r_inverse
    
    # y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] = \
        llp1 * (-density * gamma + 2. * shear_modulus * z_k2_pos * r2_inverse)
    starting_conditions_ptr[neg_index * num_ys + 1] = \
        density * ((dynamic_term / gamma) * (dynamic_term + 4. * gamma) - llp1 * gamma)
    starting_conditions_ptr[2 * num_ys + 1] = \
        2. * shear_modulus * degree_l_dbl * lm1 * r2_inverse
    
    # y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] = \
        z_k2_pos * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 2] = \
        0.
    starting_conditions_ptr[2 * num_ys + 2] = \
        r_inverse
    
    # y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] = \
        shear_modulus * (dynamic_term / beta2 - 2. * r2_inverse * z_k2_pos)
    starting_conditions_ptr[neg_index * num_ys + 3] = \
        0.
    starting_conditions_ptr[2 * num_ys + 3] = \
        2. * shear_modulus * lm1 * r2_inverse
    
    # y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] = \
        lp1 * (degree_l_dbl * gamma - dynamic_term)
    starting_conditions_ptr[neg_index * num_ys + 4] = \
        (h_k2_neg - 3.) * dynamic_term - h_k2_neg * degree_l_dbl * gamma
    starting_conditions_ptr[2 * num_ys + 4] = \
        degree_l_dbl * gamma - dynamic_term
    
    # y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[pos_index * num_ys + 4] * r_inverse
    starting_conditions_ptr[neg_index * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[neg_index * num_ys + 4] * r_inverse
    starting_conditions_ptr[2 * num_ys + 5] = \
        dlp1 * starting_conditions_ptr[2 * num_ys + 4] * r_inverse - (3. * degree_l_dbl * gamma * r_inverse)


########################################################################################################################
#### Liquid Layers
########################################################################################################################


cdef void cf_kamata_liquid_dynamic_compressible(
        double frequency,
        double radius,
        double density,
        double complex bulk_modulus,
        int degree_l,
        double G_to_use,
        size_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """  Calculate the starting guess at the bottom of a liquid layer using the dynamic assumption.

    This function uses the Kamata et al (2015; JGR:P) equations (Eq. B29-B37).

    Using the dynamic assumption in a liquid layer results in two independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), compressibility, and
       bulk dissipation (no shear dissipation within liquid layers).

    References
    ----------
    KMN15 Eq. B29-B37

    Parameters
    ----------
    frequency : double
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    radius : double
        Radius where the radial functions are calculated. [m]
    density : double
        Density at `radius` [kg m-3]
    bulk_modulus : double complex
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    degree_l : unsigned char
        Tidal harmonic order.
    G_to_use : double
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    num_ys : ssize_t
        Number of radial solutions for this layer type.
    starting_conditions_ptr : double complex*, <Output>
        Desired starting conditions for this layer.
        Two independent liquid guesses (sn1, sn2, sn3)

    """

    # Convert compressibility parameters
    # For the liquid layer the shear modulus is zero so the 1st Lame parameter = bulk modulus
    cdef double complex lame = bulk_modulus

    # Optimizations
    cdef double dynamic_term = frequency * frequency
    cdef double r_inverse    = 1. / radius
    cdef double r2           = radius * radius
    cdef double degree_l_dbl = <double>degree_l

    # Helper functions
    cdef double gamma          = (4. * d_PI_DBL * G_to_use * density / 3.)
    cdef double f              = -dynamic_term / gamma
    cdef double h              = f - (degree_l_dbl + 1.)
    cdef double complex alpha2 = lame / density
    cdef double complex k2     = (1. / alpha2) * (dynamic_term + 4. * gamma - 
                                                  (degree_l_dbl * (degree_l_dbl + 1) * gamma**2 / dynamic_term))

    # See Eqs. B33--B36 in KMN15
    # y1, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 0] = \
        -f * r_inverse * cf_z_calc(k2 * r2, degree_l=degree_l)
    starting_conditions_ptr[1 * num_ys + 0] = \
        degree_l_dbl * r_inverse
    
    # y2, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 1] = \
        -density * (f * (dynamic_term + 4 * gamma) + degree_l_dbl * (degree_l_dbl + 1) * gamma)
    starting_conditions_ptr[1 * num_ys + 1] = \
        0.
    
    # y5, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 2] = \
        3. * gamma * f - h * (degree_l_dbl * gamma - dynamic_term)
    starting_conditions_ptr[1 * num_ys + 2] = \
        degree_l_dbl * gamma - dynamic_term
    
    # y6, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 3] = \
        (2. * degree_l_dbl + 1.) * starting_conditions_ptr[0 * num_ys + 2] * r_inverse
    starting_conditions_ptr[1 * num_ys + 3] = \
        ((2. * degree_l_dbl + 1.) * starting_conditions_ptr[1 * num_ys + 2] * r_inverse) - \
        ((3. * degree_l_dbl * gamma) * r_inverse)


cdef void cf_kamata_liquid_dynamic_incompressible(
        double frequency,
        double radius,
        double density,
        int degree_l,
        double G_to_use,
        size_t num_ys, 
        double complex* starting_conditions_ptr
        ) noexcept nogil:
    """  Calculate the starting guess at the bottom of a liquid layer using the dynamic and incompressible assumption.

    This function uses the Kamata et al. (2015; JGR:P) equations (Eq. B29-B37).

    Using the dynamic assumption in a liquid layer results in two independent solutions for the radial derivatives.

    These independent solutions allow for a general tidal harmonic l, for dynamic tides (w != 0), incompressible, and
       bulk dissipation (no shear dissipation within liquid layers).

    References
    ----------
    KMN15 Eq. B29-B37

    Parameters
    ----------
    frequency : double
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    radius : double
        Radius where the radial functions are calculated. [m]
    density : double
        Density at `radius` [kg m-3]
    degree_l : unsigned char
        Tidal harmonic order.
    G_to_use : double
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    num_ys : ssize_t
        Number of radial solutions for this layer type.
    starting_conditions_ptr : double complex*, <Output>
        Desired starting conditions for this layer.
        Two independent liquid guesses (sn1, sn2, sn3)
    """

    # Optimizations
    cdef double dynamic_term = frequency * frequency
    cdef double r_inverse    = 1. / radius
    cdef double degree_l_dbl = <double>degree_l
    cdef double lp1          = degree_l_dbl + 1.
    cdef double llp1         = degree_l_dbl * lp1
    cdef double dlp1         = 2. * degree_l_dbl + 1.

    # Helper functions
    cdef double gamma = (4. * d_PI_DBL * G_to_use * density / 3.)
    cdef double f     = -dynamic_term / gamma
    cdef double h     = f - lp1

    # See Eqs. B33--B36 in KMN15
    
    # y1, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 0] = \
        0.
    starting_conditions_ptr[1 * num_ys + 0] = \
        degree_l_dbl * r_inverse
    
    # y2, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 1] = \
        -density * (f * (dynamic_term + 4 * gamma) + llp1 * gamma)
    starting_conditions_ptr[1 * num_ys + 1] = \
        0.
    
    # y5, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 2] = \
        3. * gamma * f - h * (degree_l_dbl * gamma - dynamic_term)
    starting_conditions_ptr[1 * num_ys + 2] = \
        degree_l_dbl * gamma - dynamic_term
    
    # y6, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 3] = \
        dlp1 * starting_conditions_ptr[0 * num_ys + 2] * r_inverse
    starting_conditions_ptr[1 * num_ys + 3] = \
        (dlp1 * starting_conditions_ptr[1 * num_ys + 2] * r_inverse) - ((3. * degree_l_dbl * gamma) * r_inverse)
