from typing import TYPE_CHECKING

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray

    from . import TidalPotentialModeOutput


@njit(cacheable=True)
def tidal_potential(
    radius: 'FloatArray', longitude: 'FloatArray', colatitude: 'FloatArray', time: 'FloatArray',
    orbital_frequency: 'FloatArray', eccentricity: 'FloatArray',
    host_mass: float, semi_major_axis: 'FloatArray',
    ) -> 'TidalPotentialModeOutput':
    """ Tidal gravitational potential assuming low eccentricity, no obliquity, and synchronous rotation

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
    eccentricity : FloatArray
        Eccentricity of the orbit
    host_mass : float
        Mass of tide-raising world
    semi_major_axis : FloatArray
        Orbital semi-major axis [m]

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
    """

    # Convert 'frequency' into real frequency and mode
    mode = orbital_frequency
    orbital_frequency = np.abs(orbital_frequency)

    # Associated Legendre Functions and their derivatives
    p_20 = (1. / 2.) * (3. * np.cos(colatitude)**2 - 1.)
    dp_20_dtheta = -3. * np.cos(colatitude) * np.sin(colatitude)
    dp2_20_dtheta2 = 3. * (np.sin(colatitude)**2 - np.cos(colatitude)**2)

    p_22 = 3. * (1. - np.cos(colatitude)**2)
    dp_22_dtheta = 6. * np.cos(colatitude) * np.sin(colatitude)
    dp2_22_dtheta2 = 6. * (-np.sin(colatitude)**2 + np.cos(colatitude)**2)

    # Calculate tidal potential
    # The below coefficient is often written as r2 n2 e in the literature. This assumes M_H << m. TidalPy does not
    #    use that assumption. So the n2 is replaced by G M / a3
    coeff = G * host_mass * radius**2 * eccentricity / semi_major_axis**3

    potential = coeff * \
                ((-3. / 2.) * p_20 * np.cos(orbital_frequency * time) +
                 (1. / 4.) * p_22 *
                 (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
                  4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)))

    # Its partial derivatives
    potential_partial_theta = coeff * \
                              ((-3. / 2.) * dp_20_dtheta * np.cos(orbital_frequency * time) +
                               (1. / 4.) * dp_22_dtheta *
                               (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
                                4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)))

    potential_partial_phi = coeff * \
                            ((1. / 4.) * p_22 *
                             (-6. * np.cos(orbital_frequency * time) * np.sin(2. * longitude) +
                              8. * np.sin(orbital_frequency * time) * np.cos(2. * longitude)))

    # And its 2nd order partial derivatives
    potential_partial2_theta2 = coeff * \
                                ((-3. / 2.) * dp2_20_dtheta2 * np.cos(orbital_frequency * time) +
                                 (1. / 4.) * dp2_22_dtheta2 *
                                 (3. * np.cos(orbital_frequency * time) * np.cos(2. * longitude) +
                                  4. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)))

    potential_partial2_phi2 = coeff * \
                              ((1. / 4.) * p_22 *
                               (-12. * np.cos(orbital_frequency * time) * np.cos(
                                   2. * longitude
                                   ) +  # I believe there is an error in Henning code where this (-12) is a (+12) which I believe is wrong.
                                -16. * np.sin(orbital_frequency * time) * np.sin(2. * longitude)))

    potential_partial2_theta_phi = coeff * \
                                   ((1. / 4.) * dp_22_dtheta *
                                    (-6. * np.cos(orbital_frequency * time) * np.sin(2. * longitude) +
                                     8. * np.sin(orbital_frequency * time) * np.cos(2. * longitude)))

    # Store results for this "mode"
    frequencies_by_name = {'n': orbital_frequency}
    modes_by_name = {'n': mode}
    potential_tuple_by_mode = {
        'n':
            (potential, potential_partial_theta, potential_partial_phi, potential_partial2_theta2,
             potential_partial2_phi2, potential_partial2_theta_phi)
        }

    return frequencies_by_name, modes_by_name, potential_tuple_by_mode
