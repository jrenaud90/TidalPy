# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

########################################################################################################################
#### Liquid Layers
########################################################################################################################


cdef void cf_saito_liquid_static_inccompressible(
        double radius,
        unsigned char degree_l,
        ssize_t num_ys, 
        double complex* initial_conditions_ptr
        ) noexcept nogil:
    """ Calculate the initial guess at the bottom of a liquid layer using the static assumption.

    This function uses the Saito 1974 equations (Eq. 19).

    Using the static assumption in a liquid layer results in one independent solutions for the radial derivative.

    These independent solution allow for a general tidal harmonic l, for static tides (w = 0).
    However, compressibility and all dissipation dependence is lost due to no dependence on bulk or shear moduli.


    References
    ----------
    S74

    Parameters
    ----------
    radius : float
        Radius where the radial functions are calculated. [m]
    degree_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_guesses : np.ndarray
        The one independent liquid guess (sn1)

    """

    # See Eq. 19 in Saito 1974
    # # y5 solution 0
    initial_conditions_ptr[0] = radius**degree_l

    # # y7 solution 0
    initial_conditions_ptr[1] = 2. * (degree_l - 1.) * radius**(degree_l - 1.)
