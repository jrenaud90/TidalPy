from typing import TYPE_CHECKING

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import bool_, njit

from . import MIN_SPIN_ORBITAL_DIFF

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray

    from . import TidalPotentialModeOutput


@njit(cacheable=True)
def tidal_potential(
    radius: 'FloatArray', longitude: 'FloatArray', colatitude: 'FloatArray', time: 'FloatArray',
    orbital_frequency: 'FloatArray', rotation_frequency: 'FloatArray',
    eccentricity: 'FloatArray',
    host_mass: float, semi_major_axis: 'FloatArray',
    use_static: bool = False
    ) -> 'TidalPotentialModeOutput':
    """ Tidal gravitational potential assuming moderate eccentricity, no obliquity, and non-synchronous rotation

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
    TidalPy.tides.potential.nsr_med_eccen_med_obliquity.py
    """

    # Optimizations
    cos_lat = np.cos(colatitude)
    sin_lat = np.sin(colatitude)
    cos2_lat = cos_lat * cos_lat
    sin2_lat = sin_lat * sin_lat
    cos_dbl_long = np.cos(2.0 * longitude)
    sin_dbl_long = np.sin(2.0 * longitude)

    e = eccentricity
    e2 = eccentricity * eccentricity
    e3 = eccentricity * e2
    n = orbital_frequency
    o = rotation_frequency

    # Associated Legendre Functions and their partial derivatives
    p_20 = (1. / 2.) * (3. * cos2_lat - 1.)
    dp_20_dtheta = -3. * cos_lat * sin_lat
    dp2_20_dtheta2 = 3. * (sin2_lat - cos2_lat)

    p_22 = 3. * (1. - cos2_lat)
    dp_22_dtheta = 6. * cos_lat * sin_lat
    dp2_22_dtheta2 = 6. * (cos2_lat - sin2_lat)

    zero_colat = 0. * p_20
    legendre_coeffs = (
        # P_20
        (p_20, dp_20_dtheta, dp2_20_dtheta2),
        # P_21: P_21 is unused in this truncation. set its values to zero.
        (zero_colat, zero_colat, zero_colat),
        # P_22
        (p_22, dp_22_dtheta, dp2_22_dtheta2)
        )

    # Longitude coefficients and their partial derivatives
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
        (one_long, zero_long,
         zero_long, zero_long,
         zero_long, zero_long),
        # 2*phi
        (cosine_2long_coeff, sine_2long_coeff,
         cosine_2long_coeff_dphi, sine_2long_coeff_dphi,
         cosine_2long_coeff_dphi2, sine_2long_coeff_dphi2),
        )

    # Setup mode list used under this functions assumptions.
    modes = (
        n,
        2. * n,
        3. * n,
        2. * o + n,
        2. * o - n,
        2. * o - 2. * n,
        2. * o - 3. * n,
        2. * o - 4. * n,
        2. * o - 5. * n
        )
    num_modes = 9

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
        # 2o - n
        2,
        # 2o - 2n
        2,
        # 2o - 3n
        2,
        # 2o - 4n
        2,
        # 2o - 5n
        2
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
        # 2o - n
        2,
        # 2o - 2n
        2,
        # 2o - 3n
        2,
        # 2o - 4n
        2,
        # 2o - 5n
        2
        )

    # Coefficients for each mode
    mode_coeffs = (
        # n
        -e - (9. / 8.) * e3,
        # 2n
        -(3. / 2.) * e2,
        # 3n
        -(53. / 24.) * e3,
        # 2o + n
        (1. / 288.) * e3,
        # 2o - n
        (-1. / 12.) * e + (1. / 96.) * e3,
        # 2o - 2n
        (1. / 6.) - (5. / 12.) * e2,
        # 2o - 3n
        (7. / 12.) * e - (41. / 32.) * e3,
        # 2o - 4n
        (17. / 12.) * e2,
        # 2o - 5n
        (845. / 288.) * e3
        )

    # Build storage for the potential and its derivatives
    shape = colatitude + o + e + n
    potential = np.zeros_like(shape, dtype=np.float64)
    potential_partial_theta = np.zeros_like(shape, dtype=np.float64)
    potential_partial_phi = np.zeros_like(shape, dtype=np.float64)
    potential_partial2_theta2 = np.zeros_like(shape, dtype=np.float64)
    potential_partial2_phi2 = np.zeros_like(shape, dtype=np.float64)
    potential_partial2_theta_phi = np.zeros_like(shape, dtype=np.float64)

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
        mode_switch = np.ones_like(shape, dtype=bool_)
        if not use_static:
            # Use static is False (default). Switches depend on the value of n and o
            # The orbital motion only nodes will always be on (unless n = 0 but that is not really possible).
            mode_switch *= freq > MIN_SPIN_ORBITAL_DIFF

        # Combine the mode switch with this mode's coefficient
        switch_coeff = mode_switch * mode_coeffs[mode_i]

        # Solve for the potential and its partial derivatives
        potential += switch_coeff * \
                     longitude_coeff * \
                     legendre

        potential_partial_theta += switch_coeff * \
                                   longitude_coeff * \
                                   legendre_dtheta

        potential_partial_phi += switch_coeff * \
                                 longitude_coeff_dphi * \
                                 legendre

        potential_partial2_theta2 += switch_coeff * \
                                     longitude_coeff * \
                                     legendre_dtheta2

        potential_partial2_phi2 += switch_coeff * \
                                   longitude_coeff_dphi2 * \
                                   legendre

        potential_partial2_theta_phi += switch_coeff * \
                                        longitude_coeff_dphi * \
                                        legendre_dtheta

    # # Deal with static portion of the potential
    # TODO: Is this used? It is absent from other authors definitions. For now I am including it for this function
    if use_static:
        static_coeff = (-1. / 3.) - (1. / 2.) * e2
        # The static portion for these assumptions does not depend on longitude so the partial with respect to phi is 0
        #    don't bother adding anything to those partial derivatives.
        potential += static_coeff * p_20
        potential_partial_theta += static_coeff * dp_20_dtheta
        potential_partial2_theta2 += static_coeff * dp2_20_dtheta2

    # Multiply by the outer coefficients
    global_coefficient = (3. / 2.) * G * host_mass * radius**2 / semi_major_axis**3

    potential *= global_coefficient
    potential_partial_theta *= global_coefficient
    potential_partial_phi *= global_coefficient
    potential_partial2_theta2 *= global_coefficient
    potential_partial2_phi2 *= global_coefficient
    potential_partial2_theta_phi *= global_coefficient

    # Store results for this "mode"
    # Convert 'frequency' into real frequency and mode
    mode = orbital_frequency
    orbital_frequency = np.abs(orbital_frequency)
    frequencies_by_name = {'n': orbital_frequency}
    modes_by_name = {'n': mode}
    potential_tuple_by_mode = {
        'n':
            (potential, potential_partial_theta, potential_partial_phi, potential_partial2_theta2,
             potential_partial2_phi2, potential_partial2_theta_phi)
        }

    return frequencies_by_name, modes_by_name, potential_tuple_by_mode
