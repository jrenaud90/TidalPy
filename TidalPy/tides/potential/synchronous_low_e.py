"""

"""

import numpy as np

from ...utilities.types import FloatArray

def tidal_potential(radius: float, longitude: FloatArray, colatitude: FloatArray,
                    orbital_frequency: FloatArray, eccentricity: FloatArray, time: FloatArray):

    # Associated Legendre Functions
    p_02 = (1. / 2.) * (3. * np.cos(colatitude)**2 - 1.)
    p_22 = 3. * np.sin(colatitude)**2
    dp_02_dtheta = -3. * np.cos(colatitude) * np.sin(colatitude)
    dp_22_dtheta = 6. * np.cos(colatitude) * np.sin(colatitude)

    # Calculate tidal potential
    potential = radius**2 * orbital_frequency**2 * eccentricity * \
        ( (-3. / 2.) * p_02 * np.cos(orbital_frequency * time) +
          (1. / 4.) * p_22 *
          (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
           4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)) )

    # And its partial derivatives
    potential_partial_theta = radius**2 * orbital_frequency**2 * eccentricity * \
        ( (-3. / 2.) * dp_02_dtheta * np.cos(orbital_frequency * time) +
          (1. / 4.) * dp_22_dtheta *
          (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
           4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)) )

    potential_partial_phi = radius**2 * orbital_frequency**2 * eccentricity * \
        ((-3. / 2.) * p_02 * np.cos(orbital_frequency * time) +
         (1. / 4.) * p_22 *
         (-6. * np.cos(orbital_frequency * time) * np.sin(2. * longitude) +
          8. * np.sin(orbital_frequency * time) * np.cos(2. * longitude)))
    potential_partial_phi = potential_partial_phi / np.sin(colatitude)

    return potential, potential_partial_theta, potential_partial_phi