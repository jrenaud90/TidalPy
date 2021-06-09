"""

"""

import numpy as np

from ...utilities.types import FloatArray
from ...utilities.performance import njit

@njit(cacheable=True)
def tidal_potential(radius: FloatArray, longitude: FloatArray, colatitude: FloatArray,
                    orbital_frequency: FloatArray, eccentricity: FloatArray, time: FloatArray):
    """ Tidal gravitational potential assuming low eccentricity, no obliquity, and synchronous rotation

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

    # Associated Legendre Functions and their derivatives
    p_02 = (1. / 2.) * (3. * np.cos(colatitude)**2 - 1.)
    dp_02_dtheta = -3. * np.cos(colatitude) * np.sin(colatitude)
    dp2_02_dtheta2 = 3. * (np.sin(colatitude)**2 - np.cos(colatitude)**2)

    p_22 = 3. * (1. - np.cos(colatitude)**2)
    dp_22_dtheta = 6. * np.cos(colatitude) * np.sin(colatitude)
    dp2_22_dtheta2 = 6. * (-np.sin(colatitude)**2 + np.cos(colatitude)**2)

    # Calculate tidal potential
    r2n2e = radius**2 * orbital_frequency**2 * eccentricity
    potential = r2n2e * \
        ( (-3. / 2.) * p_02 * np.cos(orbital_frequency * time) +
          (1. / 4.) * p_22 *
          (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
           4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)) )

    # Its partial derivatives
    potential_partial_theta = r2n2e  * \
        ( (-3. / 2.) * dp_02_dtheta * np.cos(orbital_frequency * time) +
          (1. / 4.) * dp_22_dtheta *
          (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
           4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)) )

    potential_partial_phi = r2n2e * \
        ( (1. / 4.) * p_22 *
          (-6. * np.cos(orbital_frequency * time) * np.sin(2. * longitude) +
           8. * np.sin(orbital_frequency * time) * np.cos(2. * longitude)))

    # And its 2nd order partial derivatives
    potential_partial2_theta2 = r2n2e * \
        ( (-3. / 2.) * dp2_02_dtheta2 * np.cos(orbital_frequency * time) +
          (1. / 4.) * dp2_22_dtheta2 *
          (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
           4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)))

    potential_partial2_phi2 = r2n2e * \
        ( (1. / 4.) * p_22 *
          (-12. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +    # I believe there is an error in Henning code where this (-12) is a (+12) which I believe is wrong.
           -16. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)))

    potential_partial2_theta_phi = r2n2e * \
        ( (1. / 4.) * dp_22_dtheta *
          (-6. * np.cos(orbital_frequency * time) * np.sin(2. * longitude) +
           8. * np.sin(orbital_frequency * time) * np.cos(2. * longitude)))



    return potential, potential_partial_theta, potential_partial_phi,\
           potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi