from typing import TYPE_CHECKING

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import bool_, njit

from . import MIN_SPIN_ORBITAL_DIFF

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray

    from . import TidalPotentialModeOutput, PotentialTupleModeOutput


@njit(cacheable=True)
def tidal_potential(
    radius: 'FloatArray', longitude: 'FloatArray', colatitude: 'FloatArray', time: 'FloatArray',
    orbital_frequency: 'FloatArray', rotation_frequency: 'FloatArray',
    eccentricity: 'FloatArray', obliquity: 'FloatArray',
    host_mass: float, semi_major_axis: 'FloatArray',
    use_static: bool = False
    ) -> 'TidalPotentialModeOutput':
    """ Tidal gravitational potential assuming moderate eccentricity, any obliquity, and non-synchronous rotation

    This function will keep the subcomponents of the potential separate for each tidal mode so that they may
       be collapsed later on (e.g., after combining with different frequency response at each mode).

    This version, which allows for obliquity, is about 2.5 times slower than the version that does not.
        TidalPy.tides.potential.nsr_modes_med_eccen_no_obliquity.py

    Parameters
    ----------
    radius : FloatArray
        Radius of the world [m]
    longitude : FloatArray
        Longitude [radians]
    colatitude : FloatArray
        Co-latitude [radians]
    time : FloatArray
        Time (orbit position) [s]
    orbital_frequency : FloatArray
        Orbital mean motion of the world [rads s-1]
    rotation_frequency: FloatArray
        Rotation rate of the planet [rad s-1]
    eccentricity : FloatArray
        Eccentricity of the orbit
    obliquity : FloatArray
        Obliquity of the planet relative to the orbital plane [rads]
    host_mass : float
        Mass of tide-raising world
    semi_major_axis : FloatArray
        Orbital semi-major axis [m]
    use_static : bool = False
        Use the static portion of the potential equation (no time dependence in these terms)

    Returns
    -------
    tidal_frequencies : Dict[str, FloatArray]
        Tidal frequencies, abs(modes) [radians s-1]
    tidal_modes : Dict[str, FloatArray]
        Tidal frequency modes [radians s-1]
    potential_tuple_by_mode : PotentialTupleModeOutput
        Storage for the tidal potential and its derivatives.

        potential : FloatArray
            Tidal Potential
        potential_partial_theta : FloatArray
            Partial Derivative of the Tidal Potential wrt colatitude.
        potential_partial_phi : FloatArray
            Partial Derivative of the Tidal Potential wrt longitude.
        potential_partial2_theta2 : FloatArray
            2nd Partial Derivative of the Tidal Potential wrt colatitude.
        potential_partial2_phi2 : FloatArray
            2nd Partial Derivative of the Tidal Potential wrt longitude.
        potential_partial2_theta_phi : FloatArray
            2nd Partial Derivative of the Tidal Potential wrt colatitude and longitude.

    See Also
    --------
    TidalPy.tides.potential.nsr_modes_med_eccen_no_obliquity.py
    """

    # Optimizations
    cos_lat = np.cos(colatitude)
    sin_lat = np.sin(colatitude)
    cos2_lat = cos_lat * cos_lat
    sin2_lat = sin_lat * sin_lat
    sqrt_1minuscos2 = np.sqrt(1. - cos2_lat)
    cos_long = np.cos(longitude)
    sin_long = np.sin(longitude)
    cos_dbl_long = np.cos(2.0 * longitude)
    sin_dbl_long = np.sin(2.0 * longitude)

    e = eccentricity
    n = orbital_frequency
    o = rotation_frequency
    obliquity_half = obliquity / 2.

    sin_o = np.sin(obliquity_half)
    sin_o_2 = sin_o * sin_o
    sin_o_3 = sin_o * sin_o_2
    sin_o_4 = sin_o * sin_o_3

    cos_o = np.cos(obliquity_half)
    cos_o_2 = cos_o * cos_o
    cos_o_3 = cos_o * cos_o_2
    cos_o_4 = cos_o * cos_o_3

    cos4_p_sin4 = cos_o_4 + sin_o_4
    cos2_sin2 = cos_o_2 * sin_o_2
    sin3_cos = sin_o_3 * cos_o
    cos3_sin = cos_o_3 * sin_o

    # Associated Legendre Functions and their partial derivatives
    p_20 = (1. / 2.) * (3. * cos2_lat - 1.)
    dp_20_dtheta = -3. * cos_lat * sin_lat
    dp2_20_dtheta2 = 3. * (sin2_lat - cos2_lat)

    p_21 = -3. * cos_lat * sqrt_1minuscos2
    dp_21_dtheta = 3. * sin_lat * (1. - 2. * cos2_lat) / sqrt_1minuscos2
    dp2_21_dtheta2 = 3. * cos_lat * (3. * sin2_lat + 2. * cos2_lat**2 - (2. * sin2_lat + 3.) * cos2_lat + 1.) / \
                     ((1. - cos2_lat) * sqrt_1minuscos2)

    p_22 = 3. * (1. - cos2_lat)
    dp_22_dtheta = 6. * cos_lat * sin_lat
    dp2_22_dtheta2 = 6. * (cos2_lat - sin2_lat)

    legendre_coeffs = (
        # P_20
        (p_20, dp_20_dtheta, dp2_20_dtheta2),
        # P_21
        (p_21, dp_21_dtheta, dp2_21_dtheta2),
        # P_22
        (p_22, dp_22_dtheta, dp2_22_dtheta2)
        )

    # Longitude coefficients and their partial derivatives
    # Note that the 1phi terms correspond to sin(phi + mode) whereas the 2phi terms are cos(2phi + mode)
    # 1\phi
    cosine_1long_coeff = sin_long
    sine_1long_coeff = cos_long
    cosine_1long_coeff_dphi = cos_long
    sine_1long_coeff_dphi = -1. * sin_long
    cosine_1long_coeff_dphi2 = -1. * sin_long
    sine_1long_coeff_dphi2 = -1. * cos_long

    # 2\phi
    cosine_2long_coeff = cos_dbl_long
    sine_2long_coeff = -sin_dbl_long
    cosine_2long_coeff_dphi = -2. * sin_dbl_long
    sine_2long_coeff_dphi = -2. * cos_dbl_long
    cosine_2long_coeff_dphi2 = -4. * cos_dbl_long
    sine_2long_coeff_dphi2 = 4. * sin_dbl_long

    zero_long = 0. * cos_dbl_long
    one_long = 1. + zero_long
    longitude_coeffs = (
        # 0*phi
        (one_long, zero_long,  # The 1 is needed in the first position otherwise the frequency dependence is lost.
         zero_long, zero_long,  # all the partials will remain zero.
         zero_long, zero_long),
        # 1*phi is unused in this truncation, set it equal to zero
        (cosine_1long_coeff, sine_1long_coeff,
         cosine_1long_coeff_dphi, sine_1long_coeff_dphi,
         cosine_1long_coeff_dphi2, sine_1long_coeff_dphi2),
        # 2*phi
        (cosine_2long_coeff, sine_2long_coeff,
         cosine_2long_coeff_dphi, sine_2long_coeff_dphi,
         cosine_2long_coeff_dphi2, sine_2long_coeff_dphi2),
        )

    # Setup mode list used under this functions assumptions.
    mode_names = (
        'n',
        '2n',
        '3n',
        '2o+n',
        '2o+2n',
        '2o+3n',
        '2o-n',
        '2o-2n',
        '2o-3n',
        'o',
        '2o',
        'o+n',
        'o+2n',
        'o+3n',
        'o-n',
        'o-2n',
        'o-3n'
        )

    # Setup mode list used under this functions assumptions.
    modes = (
        n,
        2. * n,
        3. * n,
        2. * o + n,
        2. * o + 2. * n,
        2. * o + 3. * n,
        2. * o - n,
        2. * o - 2. * n,
        2. * o - 3. * n,
        o,
        2. * o,
        o + n,
        o + 2. * n,
        o + 3. * n,
        o - n,
        o - 2. * n,
        o - 3. * n
        )
    num_modes = 17

    # Indicate which legendre polynomial is associated with which mode. The number here refers to the m integer
    mode_legendre = (
        # n
        0,
        # 2n
        0,
        # 3n
        0,
        # 2o + n
        2,
        # 2o + 2n
        2,
        # 2o + 3n
        2,
        # 2o - n
        2,
        # 2o - 2n
        2,
        # 2o - 3n
        2,
        # o
        1,
        # 2o
        2,
        # o + n
        1,
        # o + 2n
        1,
        # o + 3n
        1,
        # o - n
        1,
        # o - 2n
        1,
        # o - 3n
        1
        )

    mode_longitude = (
        # n
        0,
        # 2n
        0,
        # 3n
        0,
        # 2o + n
        2,
        # 2o + 2n
        2,
        # 2o + 3n
        2,
        # 2o - n
        2,
        # 2o - 2n
        2,
        # 2o - 3n
        2,
        # o
        1,
        # 2o
        2,
        # o + n
        1,
        # o + 2n
        1,
        # o + 3n
        1,
        # o - n
        1,
        # o - 2n
        1,
        # o - 3n
        1
        )

    # Coefficients for each mode
    mode_coeffs = (
        # n
        (-e) * cos4_p_sin4 + (5. * e) * cos2_sin2,
        # 2n
        (-2. + 11.) * cos2_sin2,
        # 3n
        (-7. * e) * cos2_sin2,
        # 2o + n
        ((-1. / 12.) * e) * sin_o_4 + ((1. / 2.) * e) * cos2_sin2,
        # 2o + 2n
        (1. / 6.) * sin_o_4,
        # 2o + 3n
        ((7. / 12.) * e) * sin_o_4,
        # 2o - n
        ((-1. / 12.) * e) * cos_o_4 + ((1. / 2.) * e) * cos2_sin2,
        # 2o - 2n
        (1. / 6.) * cos_o_4,
        # 2o - 3n
        ((7. / 12.) * e) * cos_o_4,
        # o
        (-2. / 3.) * sin3_cos,
        # 2o
        (1. / 3.) * cos2_sin2,
        # o + n
        ((-4. / 3.) * e) * sin3_cos + e * cos3_sin,
        # o + 2n
        (2. / 3.) * sin3_cos,
        # o + 3n
        ((7. / 3.) * e) * sin3_cos,
        # o - n
        (-e) * sin3_cos + ((4. / 3.) * e) * cos3_sin,
        # o - 2n
        (-2. / 3.) * cos3_sin,
        # o - 3n
        ((-7. / 3.) * e) * cos3_sin
        )

    # Prepare static coeff
    static_coeff = (-1. / 3.) * cos4_p_sin4 + (4. / 3.) * cos2_sin2
    static_term = static_coeff * p_20
    static_term_partial_theta = static_coeff * dp_20_dtheta
    static_term_partial2_theta2 = static_coeff * dp2_20_dtheta2

    # Prepare global coefficient
    global_coefficient = (3. / 2.) * G * host_mass * radius**2 / semi_major_axis**3

    # Build storage for the potential and its derivatives by mode
    oen_shape = o + e + n
    frequencies_by_name = dict()  # type: Dict[str, FloatArray]
    modes_by_name = dict()  # type: Dict[str, FloatArray]
    potential_tuple_by_mode = dict()  # type: PotentialTupleModeOutput

    # Go through modes and add their contribution to the potential and its derivatives.
    for mode_i in range(num_modes):
        # Optimizations
        mode = modes[mode_i]
        cos_mode = np.cos(mode * time)
        sin_mode = np.sin(mode * time)
        freq = np.abs(mode)

        # Pull out legendre coeffs for this mode
        legendre, legendre_dtheta, legendre_dtheta2 = legendre_coeffs[mode_legendre[mode_i]]

        # Pull out longitude coeffs for this mode
        cosine_coeff, sine_coeff, \
            cosine_coeff_dphi, sine_coeff_dphi, \
            cosine_coeff_dphi2, sine_coeff_dphi2 = longitude_coeffs[mode_longitude[mode_i]]
        longitude_coeff = (cosine_coeff * cos_mode + sine_coeff * sin_mode)
        longitude_coeff_dphi = (cosine_coeff_dphi * cos_mode + sine_coeff_dphi * sin_mode)
        longitude_coeff_dphi2 = (cosine_coeff_dphi2 * cos_mode + sine_coeff_dphi2 * sin_mode)

        # Switches
        # There will be terms that are non-zero even though they do not carry a time dependence. This switch will
        #   ensure all non-time dependence --> zero unless the user sets `use_static` = True.
        mode_switch = np.ones_like(oen_shape, dtype=bool_)
        if not use_static:
            # Use static is False (default). Switches depend on the value of n and o
            # The orbital motion only nodes will always be on (unless n = 0 but that is not really possible).
            mode_switch *= freq > MIN_SPIN_ORBITAL_DIFF

        # Combine the mode switch with this mode's coefficient
        switch_coeff = mode_switch * mode_coeffs[mode_i]

        # Solve for the potential and its partial derivatives
        potential = switch_coeff * \
                    longitude_coeff * \
                    legendre

        potential_partial_theta = switch_coeff * \
                                  longitude_coeff * \
                                  legendre_dtheta

        potential_partial_phi = switch_coeff * \
                                longitude_coeff_dphi * \
                                legendre

        potential_partial2_theta2 = switch_coeff * \
                                    longitude_coeff * \
                                    legendre_dtheta2

        potential_partial2_phi2 = switch_coeff * \
                                  longitude_coeff_dphi2 * \
                                  legendre

        potential_partial2_theta_phi = switch_coeff * \
                                       longitude_coeff_dphi * \
                                       legendre_dtheta

        # # Deal with static portion of the potential
        # TODO: Is this used? It is absent from other authors definitions. For now I am including it for this function
        if use_static:
            # The static portion for these assumptions does not depend on longitude so the partial with respect to phi is 0
            #    don't bother adding anything to those partial derivatives.
            potential += static_term
            potential_partial_theta += static_term_partial_theta
            potential_partial2_theta2 += static_term_partial2_theta2

        # Multiply by the outer coefficients
        potential *= global_coefficient
        potential_partial_theta *= global_coefficient
        potential_partial_phi *= global_coefficient
        potential_partial2_theta2 *= global_coefficient
        potential_partial2_phi2 *= global_coefficient
        potential_partial2_theta_phi *= global_coefficient

        # Store results for this mode
        mode_name = mode_names[mode_i]
        frequencies_by_name[mode_name] = freq
        modes_by_name[mode_name] = mode
        potential_tuple_by_mode[mode_name] = \
            (potential, potential_partial_theta, potential_partial_phi, potential_partial2_theta2,
             potential_partial2_phi2, potential_partial2_theta_phi)

    return frequencies_by_name, modes_by_name, potential_tuple_by_mode
