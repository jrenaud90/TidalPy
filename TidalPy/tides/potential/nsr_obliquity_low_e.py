import numpy as np

from . import TidalPotentialOutput
from ...utilities.performance import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def tidal_potential(
    radius: FloatArray, longitude: FloatArray, colatitude: FloatArray,
    orbital_frequency: FloatArray, eccentricity: FloatArray, time: FloatArray,
    obliquity: FloatArray, rotation_rate: FloatArray, periapsis: float, world_radius: float
    ) -> TidalPotentialOutput:
    """ Tidal gravitational potential assuming low eccentricity, no obliquity, and non-synchronous rotation

    Parameters
    ----------
    radius : FloatArray
        Radius of the world [m]
    longitude : FloatArray
        Longitude [radians]
    colatitude : FloatArray
        Co-latitude [radians]
    orbital_frequency : FloatArray
        Orbital mean motion of the world [rads s-1]
    eccentricity : FloatArray
        Eccentricity of the orbit
    time : FloatArray
        Time (orbit position) [s]
    obliquity : FloatArray
        World's obliquity [radians]
    rotation_rate: FloatArray
        Rotation rate of the planet [rad s-1]
    periapsis: float
        Argument of periapsis (used with rotation rate as a phase offset) [radians]
    world_radius: float
        World's surface radius [m]

    Returns
    -------
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
    """

    # Optimizations
    cos_lat = np.cos(colatitude)
    sin_lat = np.sin(colatitude)
    cos2_lat = cos_lat * cos_lat
    sin2_lat = sin_lat * sin_lat
    sqrt_cos2_lat = np.sqrt(1. - cos2_lat)

    # Associated Legendre Functions and their derivatives
    p_20 = (1. / 2.) * (3. * cos2_lat - 1.)
    dp_20_dtheta = -3. * cos_lat * sin_lat
    dp2_20_dtheta2 = 3. * (sin2_lat - cos2_lat)

    p_21 = -3. * cos_lat * sqrt_cos2_lat
    dp_21_dtheta = (3. * sin_lat * (1. - 2. * cos2_lat)) / sqrt_cos2_lat
    dp2_21_dtheta2 = 3. * cos_lat * sqrt_cos2_lat + \
                     6. * sin2_lat * cos_lat / sqrt_cos2_lat - \
                     3. * cos_lat * ((cos2_lat / sqrt_cos2_lat) - (sin2_lat * cos2_lat / ((1. - cos2_lat)**(3 / 2))) -
                                     sin2_lat / sqrt_cos2_lat)

    p_22 = 3. * (1. - cos2_lat)
    dp_22_dtheta = 6. * cos_lat * sin_lat
    dp2_22_dtheta2 = 6. * (cos2_lat - sin2_lat)

    # Radius Factor
    # # TODO: Most authors use (R^2), but (R^5/r^3) shows up in the derivation. Should the lower case r be included?
    #       Looking at Wahr+ (2009; Icarus) it looks like r^2 is the correct choice. There is a (r/R)^2 * R^2
    # radius_factor = world_radius**2 * np.ones_like(radius)
    radius_factor = radius**2
    # radius_factor = (world_radius**5 / radius**3)

    # Calculate the sub-components of the tidal potential
    # See Eqs. 2 -- 7 in Jara-Orué & Vermeerson (2011; Icarus)

    # # Static
    # # TODO: Is this used? It is absent from other authors definitions. For now I am excluding it
    # compo_static = \
    #     -(1. / 2.) * p_20 + (1. / 4.) * p_22 * np.cos(2. * longitude)
    # compo_static_partial_theta = \
    #     -(1. / 2.) * dp_20_dtheta + (1. / 4.) * dp_22_dtheta * np.cos(2. * longitude)
    # compo_static_partial_phi = \
    #     -(1. / 2.) * p_22 * np.sin(2. * longitude)
    # compo_static_partial2_theta2 = \
    #     -(1. / 2.) * dp2_20_dtheta2 + (1. / 4.) * dp2_22_dtheta2 * np.cos(2. * longitude)
    # compo_static_partial2_phi2 = \
    #     -p_22 * np.cos(2. * longitude)
    # compo_static_partial2_theta_phi = \
    #     -(1. / 2.) * dp_22_dtheta * np.sin(2. * longitude)

    # # NSR
    nsr_factor = rotation_rate - orbital_frequency
    compo_nsr = \
        -(1. / 2.) * p_22 * np.sin(2. * longitude + nsr_factor * time) * np.sin(nsr_factor * time)
    compo_nsr_partial_theta = \
        -(1. / 2.) * dp_22_dtheta * np.sin(2. * longitude + nsr_factor * time) * np.sin(nsr_factor * time)
    compo_nsr_partial_phi = \
        -p_22 * np.cos(2. * longitude + nsr_factor * time) * np.sin(nsr_factor * time)
    compo_nsr_partial2_theta2 = \
        -(1. / 2.) * dp2_22_dtheta2 * np.sin(2. * longitude + nsr_factor * time) * np.sin(nsr_factor * time)
    compo_nsr_partial2_phi2 = \
        2. * p_22 * np.sin(2. * longitude + nsr_factor * time) * np.sin(nsr_factor * time)
    compo_nsr_partial2_theta_phi = \
        -dp_22_dtheta * np.cos(2. * longitude + nsr_factor * time) * np.sin(nsr_factor * time)

    # # Eccentricity Term 1
    compo_e1 = \
        -(3. / 2.) * eccentricity * p_20 * np.cos(orbital_frequency * time)
    compo_e1_partial_theta = \
        -(3. / 2.) * eccentricity * dp_20_dtheta * np.cos(orbital_frequency * time)
    compo_e1_partial_phi = \
        0.0
    compo_e1_partial2_theta2 = \
        -(3. / 2.) * eccentricity * dp2_20_dtheta2 * np.cos(orbital_frequency * time)
    compo_e1_partial2_phi2 = \
        0.0
    compo_e1_partial2_theta_phi = \
        0.0

    # # Eccentricity Term 2
    compo_e2 = \
        (eccentricity / 4.) * p_22 * \
        (3. * np.cos(2. * longitude) * np.cos(orbital_frequency * time) +
         4. * np.sin(2. * longitude) * np.sin(orbital_frequency * time))
    compo_e2_partial_theta = \
        (eccentricity / 4.) * dp_22_dtheta * \
        (3. * np.cos(2. * longitude) * np.cos(orbital_frequency * time) +
         4. * np.sin(2. * longitude) * np.sin(orbital_frequency * time))
    compo_e2_partial_phi = \
        (eccentricity / 4.) * p_22 * \
        (-6. * np.sin(2. * longitude) * np.cos(orbital_frequency * time) +
         8. * np.cos(2. * longitude) * np.sin(orbital_frequency * time))
    compo_e2_partial2_theta2 = \
        (eccentricity / 4.) * dp2_22_dtheta2 * \
        (3. * np.cos(2. * longitude) * np.cos(orbital_frequency * time) +
         4. * np.sin(2. * longitude) * np.sin(orbital_frequency * time))
    compo_e2_partial2_phi2 = \
        (eccentricity / 4.) * p_22 * \
        (-12. * np.cos(2. * longitude) * np.cos(orbital_frequency * time) +
         -16. * np.sin(2. * longitude) * np.sin(orbital_frequency * time))
    # I believe there is an error in Henning code where this^^ (-12) is a (+12) which I believe is wrong.
    compo_e2_partial2_theta_phi = \
        (eccentricity / 4.) * dp_22_dtheta * \
        (-6. * np.sin(2. * longitude) * np.cos(orbital_frequency * time) +
         8. * np.cos(2. * longitude) * np.sin(orbital_frequency * time))

    # # Obliquity Term
    compo_obli = \
        p_21 * np.cos(obliquity) * np.sin(obliquity) * np.cos(longitude) * \
        np.sin(periapsis + orbital_frequency * time)
    compo_obli_partial_theta = \
        dp_21_dtheta * np.cos(obliquity) * np.sin(obliquity) * np.cos(longitude) * \
        np.sin(periapsis + orbital_frequency * time)
    compo_obli_partial_phi = \
        -p_21 * np.cos(obliquity) * np.sin(obliquity) * np.sin(longitude) * \
        np.sin(periapsis + orbital_frequency * time)
    compo_obli_partial2_theta2 = \
        dp2_21_dtheta2 * np.cos(obliquity) * np.sin(obliquity) * np.cos(longitude) * \
        np.sin(periapsis + orbital_frequency * time)
    compo_obli_partial2_phi2 = \
        -p_21 * np.cos(obliquity) * np.sin(obliquity) * np.cos(longitude) * \
        np.sin(periapsis + orbital_frequency * time)
    compo_obli_partial2_theta_phi = \
        -dp_21_dtheta * np.cos(obliquity) * np.sin(obliquity) * np.sin(longitude) * \
        np.sin(periapsis + orbital_frequency * time)

    # Final Potential and Potential Derivatives
    coefficient = radius_factor * orbital_frequency**2
    potential = coefficient * (compo_nsr + compo_e1 + compo_e2 + compo_obli)

    potential_partial_theta = coefficient * \
                              (compo_nsr_partial_theta + compo_e1_partial_theta +
                               compo_e2_partial_theta + compo_obli_partial_theta)

    potential_partial_phi = coefficient * \
                            (compo_nsr_partial_phi + compo_e1_partial_phi +
                             compo_e2_partial_phi + compo_obli_partial_phi)

    potential_partial2_theta2 = coefficient * \
                                (compo_nsr_partial2_theta2 + compo_e1_partial2_theta2 +
                                 compo_e2_partial2_theta2 + compo_obli_partial2_theta2)

    potential_partial2_phi2 = coefficient * \
                              (compo_nsr_partial2_phi2 + compo_e1_partial2_phi2 +
                               compo_e2_partial2_phi2 + compo_obli_partial2_phi2)

    potential_partial2_theta_phi = coefficient * \
                                   (compo_nsr_partial2_theta_phi + compo_e1_partial2_theta_phi +
                                    compo_e2_partial2_theta_phi + compo_obli_partial2_theta_phi)

    return potential, potential_partial_theta, potential_partial_phi, \
           potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi