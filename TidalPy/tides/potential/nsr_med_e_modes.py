from typing import Dict, Tuple

import numpy as np

from . import MIN_SPIN_ORBITAL_DIFF
from ...constants import G
from ...utilities.performance import bool_, njit
from ...utilities.types import FloatArray

TidalPotentialModeOutput = Tuple[Dict[str, FloatArray], Dict[str, FloatArray], Dict[str, FloatArray],
                                 Dict[str, FloatArray], Dict[str, FloatArray], Dict[str, FloatArray],
                                 Dict[str, FloatArray]]


@njit(cacheable=True)
def tidal_potential(
    radius: FloatArray, longitude: FloatArray, colatitude: FloatArray,
    orbital_frequency: FloatArray, eccentricity: FloatArray, time: FloatArray,
    rotation_rate: FloatArray, world_radius: float, host_mass: float, semi_major_axis: FloatArray,
    use_static: bool = False,
    ) -> TidalPotentialModeOutput:
    """ Tidal gravitational potential assuming moderate eccentricity, no obliquity, and non-synchronous rotation
    utilizing multiple tidal modes.

    This method uses multiple frequency models that are linear combinations of orbital motion and planetary rotation
       rate.

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
    rotation_rate: FloatArray
        Rotation rate of the planet [rad s-1]
    world_radius: float
        World's surface radius [m]
    host_mass : float
        Mass of tide-rasing world
    semi_major_axis : FloatArray
        Orbital semi-major axis [m]
    use_static : bool = False
        Use the static portion of the potential equation (no time dependence in these terms)

    Returns
    -------
    tidal_modes : Dict[str, FloatArray]
        Tidal frequency modes [radians s-1]
    potential : Dict[str, FloatArray]
        Tidal Potential
    potential_partial_theta : Dict[str, FloatArray]
        Partial Derivative of the Tidal Potential wrt colatitude.
    potential_partial_phi : Dict[str, FloatArray]
        Partial Derivative of the Tidal Potential wrt longitude.
    potential_partial2_theta2 : Dict[str, FloatArray]
        2nd Partial Derivative of the Tidal Potential wrt colatitude.
    potential_partial2_phi2 : Dict[str, FloatArray]
        2nd Partial Derivative of the Tidal Potential wrt longitude.
    potential_partial2_theta_phi : Dict[str, FloatArray]
        2nd Partial Derivative of the Tidal Potential wrt colatitude and longitude.
    """

    # Optimizations
    cos_lat = np.cos(colatitude)
    sin_lat = np.sin(colatitude)
    cos2_lat = cos_lat * cos_lat
    sin2_lat = sin_lat * sin_lat
    sqrt_cos2_lat = np.sqrt(1. - cos2_lat)
    dbl_long = 2.0 * longitude
    e = eccentricity
    e2 = eccentricity * eccentricity
    e3 = eccentricity * e2
    n = orbital_frequency
    o = rotation_rate
    shape = o * n * time * colatitude * eccentricity

    # Associated Legendre Functions and their derivatives
    p_20 = (1. / 2.) * (3. * cos2_lat - 1.)
    dp_20_dtheta = -3. * cos_lat * sin_lat
    dp2_20_dtheta2 = 3. * (sin2_lat - cos2_lat)

    p_22 = 3. * (1. - cos2_lat)
    dp_22_dtheta = 6. * cos_lat * sin_lat
    dp2_22_dtheta2 = 6. * (cos2_lat - sin2_lat)

    # # Build frequency modes
    modes = {
        'n'    : np.abs(n),
        '2n'   : np.abs(2. * n),
        '3n'   : np.abs(3. * n),
        '2o+n' : np.abs(2. * o + n),
        '2o-n' : np.abs(2. * o - n),
        '2o-2n': np.abs(2. * o - 2. * n),
        '2o-3n': np.abs(2. * o - 3. * n),
        '2o-4n': np.abs(2. * o - 4. * n),
        '2o-5n': np.abs(2. * o - 5. * n)
        }

    # There will be terms that are non-zero even though they do not carry a time dependence. This switch will
    #   ensure all non-time dependence --> zero unless the user sets `use_static` = True.
    mode_switch = {
        'n'    : np.ones_like(shape, dtype=bool_),
        '2n'   : np.ones_like(shape, dtype=bool_),
        '3n'   : np.ones_like(shape, dtype=bool_),
        '2o+n' : np.ones_like(shape, dtype=bool_),
        '2o-n' : np.ones_like(shape, dtype=bool_),
        '2o-2n': np.ones_like(shape, dtype=bool_),
        '2o-3n': np.ones_like(shape, dtype=bool_),
        '2o-4n': np.ones_like(shape, dtype=bool_),
        '2o-5n': np.ones_like(shape, dtype=bool_)
        }
    if use_static:
        # Use static is True. All switches are on regardless of n and o
        pass
    else:
        # Use static is False (default). Switches depend on the value of n and o
        # The orbital motion only nodes will always be on (unless n = 0 but that is not really possible).
        mode_switch['2o+n'] *= np.abs(2. * o + n) > MIN_SPIN_ORBITAL_DIFF
        mode_switch['2o-n'] *= np.abs(2. * o - n) > MIN_SPIN_ORBITAL_DIFF
        mode_switch['2o-2n'] *= np.abs(2. * o - 2. * n) > MIN_SPIN_ORBITAL_DIFF
        mode_switch['2o-3n'] *= np.abs(2. * o - 3. * n) > MIN_SPIN_ORBITAL_DIFF
        mode_switch['2o-4n'] *= np.abs(2. * o - 4. * n) > MIN_SPIN_ORBITAL_DIFF
        mode_switch['2o-5n'] *= np.abs(2. * o - 5. * n) > MIN_SPIN_ORBITAL_DIFF

    # # Static
    # # TODO: Is this used? It is absent from other authors definitions. For now I am including it for this function
    if use_static:
        compo_static = \
            -(1. / 3.) * p_20 * (1. / 9.)  # The 1/9 is to spread it over active modes.
        compo_static_partial_theta = \
            -(1. / 3.) * dp_20_dtheta * (1. / 9.)  # The 1/9 is to spread it over active modes.
        compo_static_partial_phi = \
            np.zeros_like(p_20)
        compo_static_partial2_theta2 = \
            -(1. / 3.) * dp2_20_dtheta2 * (1. / 9.)  # The 1/9 is to spread it over active modes.
        compo_static_partial2_phi2 = \
            np.zeros_like(p_20)
        compo_static_partial2_theta_phi = \
            np.zeros_like(p_20)
    else:
        compo_static = \
            np.zeros_like(p_20)
        compo_static_partial_theta = \
            np.zeros_like(p_20)
        compo_static_partial_phi = \
            np.zeros_like(p_20)
        compo_static_partial2_theta2 = \
            np.zeros_like(p_20)
        compo_static_partial2_phi2 = \
            np.zeros_like(p_20)
        compo_static_partial2_theta_phi = \
            np.zeros_like(p_20)

    # # NSR Solo Term
    # Set all to zero for now.
    compo_nsr = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': (1. / 6.) * p_22 * np.cos(dbl_long + 2. * (o - n) * time) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_nsr_partial_theta = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': (1. / 6.) * dp_22_dtheta * np.cos(dbl_long + 2. * (o - n) * time) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_nsr_partial_phi = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': (-1. / 3.) * p_22 * np.sin(dbl_long + 2. * (o - n) * time) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_nsr_partial2_theta2 = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': (1. / 6.) * dp2_22_dtheta2 * np.cos(dbl_long + 2. * (o - n) * time) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_nsr_partial2_phi2 = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': (-2. / 3.) * p_22 * np.cos(dbl_long + 2. * (o - n) * time) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_nsr_partial2_theta_phi = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': (-1. / 3.) * dp_22_dtheta * np.sin(dbl_long + 2. * (o - n) * time) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }

    # # Eccentricity Solo Term
    if use_static:
        e2_static = -(1. / 2.) * (1. / 3.)  # The 1/3 is to spread it over active modes.
    else:
        e2_static = 0.0

    compo_e = {
        'n'    : p_20 * ((e2 * e2_static) + np.cos(n * time) * (-e + (-9. / 8.) * e3)) *
                 np.ones_like(shape) * mode_switch['n'],
        '2n'   : p_20 * ((e2 * e2_static) + np.cos(2. * n * time) * ((-3. / 2.) * e2)) *
                 np.ones_like(shape) * mode_switch['2n'],
        '3n'   : p_20 * ((e2 * e2_static) + np.cos(3. * n * time) * ((-53. / 24.) * e3)) *
                 np.ones_like(shape) * mode_switch['3n'],
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': np.zeros_like(shape),
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_e_partial_theta = {
        'n'    : dp_20_dtheta * ((e2 * e2_static) + np.cos(n * time) * (-e + (-9. / 8.) * e3)) *
                 np.ones_like(shape) * mode_switch['n'],
        '2n'   : dp_20_dtheta * ((e2 * e2_static) + np.cos(2. * n * time) * ((-3. / 2.) * e2)) *
                 np.ones_like(shape) * mode_switch['2n'],
        '3n'   : dp_20_dtheta * ((e2 * e2_static) + np.cos(3. * n * time) * ((-53. / 24.) * e3)) *
                 np.ones_like(shape) * mode_switch['3n'],
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': np.zeros_like(shape),
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_e_partial_phi = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': np.zeros_like(shape),
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_e_partial2_theta2 = {
        'n'    : dp2_20_dtheta2 * ((e2 * e2_static) + np.cos(n * time) * (-e + (-9. / 8.) * e3)) *
                 np.ones_like(shape) * mode_switch['n'],
        '2n'   : dp2_20_dtheta2 * ((e2 * e2_static) + np.cos(2. * n * time) * ((-3. / 2.) * e2)) *
                 np.ones_like(shape) * mode_switch['2n'],
        '3n'   : dp2_20_dtheta2 * ((e2 * e2_static) + np.cos(3. * n * time) * ((-53. / 24.) * e3)) *
                 np.ones_like(shape) * mode_switch['3n'],
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': np.zeros_like(shape),
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_e_partial2_phi2 = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': np.zeros_like(shape),
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }
    compo_e_partial2_theta_phi = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : np.zeros_like(shape),
        '2o-n' : np.zeros_like(shape),
        '2o-2n': np.zeros_like(shape),
        '2o-3n': np.zeros_like(shape),
        '2o-4n': np.zeros_like(shape),
        '2o-5n': np.zeros_like(shape)
        }

    # # Eccentricity / NSR Cross Terms
    compo_e_nsr = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : p_22 * (np.cos(dbl_long + 2. * (o + (1. / 2.) * n) * time) * ((1. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o+n'],
        '2o-n' : p_22 * (np.cos(dbl_long + 2. * (o - (1. / 2.) * n) * time) * ((-1. / 12.) * e + (1. / 96.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-n'],
        '2o-2n': p_22 * (np.cos(dbl_long + 2. * (o - n) * time) * ((-5. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': p_22 * (np.cos(dbl_long + 2. * (o - (3. / 2.) * n) * time) * ((7. / 12.) * e - (41. / 32.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-3n'],
        '2o-4n': p_22 * (np.cos(dbl_long + 2. * (o - 2. * n) * time) * ((17. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-4n'],
        '2o-5n': p_22 * (np.cos(dbl_long + 2. * (o - (5. / 2.) * n) * time) * ((845. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-5n']
        }
    compo_e_nsr_partial_theta = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : dp_22_dtheta * (np.cos(dbl_long + 2. * (o + (1. / 2.) * n) * time) * ((1. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o+n'],
        '2o-n' : dp_22_dtheta * (
                np.cos(dbl_long + 2. * (o - (1. / 2.) * n) * time) * ((-1. / 12.) * e + (1. / 96.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-n'],
        '2o-2n': dp_22_dtheta * (np.cos(dbl_long + 2. * (o - n) * time) * ((-5. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': dp_22_dtheta * (
                np.cos(dbl_long + 2. * (o - (3. / 2.) * n) * time) * ((7. / 12.) * e - (41. / 32.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-3n'],
        '2o-4n': dp_22_dtheta * (np.cos(dbl_long + 2. * (o - 2. * n) * time) * ((17. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-4n'],
        '2o-5n': dp_22_dtheta * (np.cos(dbl_long + 2. * (o - (5. / 2.) * n) * time) * ((845. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-5n']
        }
    compo_e_nsr_partial_phi = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : p_22 * (-2. * np.sin(dbl_long + 2. * (o + (1. / 2.) * n) * time) * ((1. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o+n'],
        '2o-n' : p_22 * (
                -2. * np.sin(dbl_long + 2. * (o - (1. / 2.) * n) * time) * ((-1. / 12.) * e + (1. / 96.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-n'],
        '2o-2n': p_22 * (-2. * np.sin(dbl_long + 2. * (o - n) * time) * ((-5. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': p_22 * (
                -2. * np.sin(dbl_long + 2. * (o - (3. / 2.) * n) * time) * ((7. / 12.) * e - (41. / 32.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-3n'],
        '2o-4n': p_22 * (-2. * np.sin(dbl_long + 2. * (o - 2. * n) * time) * ((17. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-4n'],
        '2o-5n': p_22 * (-2. * np.sin(dbl_long + 2. * (o - (5. / 2.) * n) * time) * ((845. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-5n']
        }
    compo_e_nsr_partial2_theta2 = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : dp2_22_dtheta2 * (np.cos(dbl_long + 2. * (o + (1. / 2.) * n) * time) * ((1. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o+n'],
        '2o-n' : dp2_22_dtheta2 * (
                np.cos(dbl_long + 2. * (o - (1. / 2.) * n) * time) * ((-1. / 12.) * e + (1. / 96.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-n'],
        '2o-2n': dp2_22_dtheta2 * (np.cos(dbl_long + 2. * (o - n) * time) * ((-5. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': dp2_22_dtheta2 * (
                np.cos(dbl_long + 2. * (o - (3. / 2.) * n) * time) * ((7. / 12.) * e - (41. / 32.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-3n'],
        '2o-4n': dp2_22_dtheta2 * (np.cos(dbl_long + 2. * (o - 2. * n) * time) * ((17. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-4n'],
        '2o-5n': dp2_22_dtheta2 * (np.cos(dbl_long + 2. * (o - (5. / 2.) * n) * time) * ((845. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-5n']
        }
    compo_e_nsr_partial2_phi2 = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : p_22 * (-4. * np.cos(dbl_long + 2. * (o + (1. / 2.) * n) * time) * ((1. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o+n'],
        '2o-n' : p_22 * (
                -4. * np.cos(dbl_long + 2. * (o - (1. / 2.) * n) * time) * ((-1. / 12.) * e + (1. / 96.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-n'],
        '2o-2n': p_22 * (-4. * np.cos(dbl_long + 2. * (o - n) * time) * ((-5. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': p_22 * (
                -4. * np.cos(dbl_long + 2. * (o - (3. / 2.) * n) * time) * ((7. / 12.) * e - (41. / 32.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-3n'],
        '2o-4n': p_22 * (-4. * np.cos(dbl_long + 2. * (o - 2. * n) * time) * ((17. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-4n'],
        '2o-5n': p_22 * (-4. * np.cos(dbl_long + 2. * (o - (5. / 2.) * n) * time) * ((845. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-5n']
        }
    compo_e_nsr_partial2_theta_phi = {
        'n'    : np.zeros_like(shape),
        '2n'   : np.zeros_like(shape),
        '3n'   : np.zeros_like(shape),
        '2o+n' : dp_22_dtheta * (-2. * np.sin(dbl_long + 2. * (o + (1. / 2.) * n) * time) * ((1. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o+n'],
        '2o-n' : dp_22_dtheta * (
                -2. * np.sin(dbl_long + 2. * (o - (1. / 2.) * n) * time) * ((-1. / 12.) * e + (1. / 96.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-n'],
        '2o-2n': dp_22_dtheta * (-2. * np.sin(dbl_long + 2. * (o - n) * time) * ((-5. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-2n'],
        '2o-3n': dp_22_dtheta * (
                -2. * np.sin(dbl_long + 2. * (o - (3. / 2.) * n) * time) * ((7. / 12.) * e - (41. / 32.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-3n'],
        '2o-4n': dp_22_dtheta * (-2. * np.sin(dbl_long + 2. * (o - 2. * n) * time) * ((17. / 12.) * e2)) *
                 np.ones_like(shape) * mode_switch['2o-4n'],
        '2o-5n': dp_22_dtheta * (-2. * np.sin(dbl_long + 2. * (o - (5. / 2.) * n) * time) * ((845. / 288.) * e3)) *
                 np.ones_like(shape) * mode_switch['2o-5n']
        }

    # Final Potential and Potential Derivatives
    susceptibility_reduced = (3. / 2.) * G * host_mass * world_radius**2 / semi_major_axis**3
    radius_factor = (radius / world_radius)**2
    coefficient = susceptibility_reduced * radius_factor

    potential = {
        'n'    : coefficient * (compo_static + compo_nsr['n'] +
                                compo_e['n'] + compo_e_nsr['n']),
        '2n'   : coefficient * (compo_static + compo_nsr['2n'] +
                                compo_e['2n'] + compo_e_nsr['2n']),
        '3n'   : coefficient * (compo_static + compo_nsr['3n'] +
                                compo_e['3n'] + compo_e_nsr['3n']),
        '2o+n' : coefficient * (compo_static + compo_nsr['2o+n'] +
                                compo_e['2o+n'] + compo_e_nsr['2o+n']),
        '2o-n' : coefficient * (compo_static + compo_nsr['2o-n'] +
                                compo_e['2o-n'] + compo_e_nsr['2o-n']),
        '2o-2n': coefficient * (compo_static + compo_nsr['2o-2n'] +
                                compo_e['2o-2n'] + compo_e_nsr['2o-2n']),
        '2o-3n': coefficient * (compo_static + compo_nsr['2o-3n'] +
                                compo_e['2o-3n'] + compo_e_nsr['2o-3n']),
        '2o-4n': coefficient * (compo_static + compo_nsr['2o-4n'] +
                                compo_e['2o-4n'] + compo_e_nsr['2o-4n']),
        '2o-5n': coefficient * (compo_static + compo_nsr['2o-5n'] +
                                compo_e['2o-5n'] + compo_e_nsr['2o-5n'])
        }

    potential_partial_theta = {
        'n'    : coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['n'] +
                                compo_e_partial_theta['n'] + compo_e_nsr_partial_theta['n']),
        '2n'   : coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2n'] +
                                compo_e_partial_theta['2n'] + compo_e_nsr_partial_theta['2n']),
        '3n'   : coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['3n'] +
                                compo_e_partial_theta['3n'] + compo_e_nsr_partial_theta['3n']),
        '2o+n' : coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2o+n'] +
                                compo_e_partial_theta['2o+n'] + compo_e_nsr_partial_theta['2o+n']),
        '2o-n' : coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2o-n'] +
                                compo_e_partial_theta['2o-n'] + compo_e_nsr_partial_theta['2o-n']),
        '2o-2n': coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2o-2n'] +
                                compo_e_partial_theta['2o-2n'] + compo_e_nsr_partial_theta['2o-2n']),
        '2o-3n': coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2o-3n'] +
                                compo_e_partial_theta['2o-3n'] + compo_e_nsr_partial_theta['2o-3n']),
        '2o-4n': coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2o-4n'] +
                                compo_e_partial_theta['2o-4n'] + compo_e_nsr_partial_theta['2o-4n']),
        '2o-5n': coefficient * (compo_static_partial_theta + compo_nsr_partial_theta['2o-5n'] +
                                compo_e_partial_theta['2o-5n'] + compo_e_nsr_partial_theta['2o-5n'])
        }

    potential_partial_phi = {
        'n'    : coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['n'] +
                                compo_e_partial_phi['n'] + compo_e_nsr_partial_phi['n']),
        '2n'   : coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2n'] +
                                compo_e_partial_phi['2n'] + compo_e_nsr_partial_phi['2n']),
        '3n'   : coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['3n'] +
                                compo_e_partial_phi['3n'] + compo_e_nsr_partial_phi['3n']),
        '2o+n' : coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2o+n'] +
                                compo_e_partial_phi['2o+n'] + compo_e_nsr_partial_phi['2o+n']),
        '2o-n' : coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2o-n'] +
                                compo_e_partial_phi['2o-n'] + compo_e_nsr_partial_phi['2o-n']),
        '2o-2n': coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2o-2n'] +
                                compo_e_partial_phi['2o-2n'] + compo_e_nsr_partial_phi['2o-2n']),
        '2o-3n': coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2o-3n'] +
                                compo_e_partial_phi['2o-3n'] + compo_e_nsr_partial_phi['2o-3n']),
        '2o-4n': coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2o-4n'] +
                                compo_e_partial_phi['2o-4n'] + compo_e_nsr_partial_phi['2o-4n']),
        '2o-5n': coefficient * (compo_static_partial_phi + compo_nsr_partial_phi['2o-5n'] +
                                compo_e_partial_phi['2o-5n'] + compo_e_nsr_partial_phi['2o-5n'])
        }

    potential_partial2_theta2 = {
        'n'    : coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['n'] +
                                compo_e_partial2_theta2['n'] + compo_e_nsr_partial2_theta2['n']),
        '2n'   : coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2n'] +
                                compo_e_partial2_theta2['2n'] + compo_e_nsr_partial2_theta2['2n']),
        '3n'   : coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['3n'] +
                                compo_e_partial2_theta2['3n'] + compo_e_nsr_partial2_theta2['3n']),
        '2o+n' : coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2o+n'] +
                                compo_e_partial2_theta2['2o+n'] + compo_e_nsr_partial2_theta2['2o+n']),
        '2o-n' : coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2o-n'] +
                                compo_e_partial2_theta2['2o-n'] + compo_e_nsr_partial2_theta2['2o-n']),
        '2o-2n': coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2o-2n'] +
                                compo_e_partial2_theta2['2o-2n'] + compo_e_nsr_partial2_theta2['2o-2n']),
        '2o-3n': coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2o-3n'] +
                                compo_e_partial2_theta2['2o-3n'] + compo_e_nsr_partial2_theta2['2o-3n']),
        '2o-4n': coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2o-4n'] +
                                compo_e_partial2_theta2['2o-4n'] + compo_e_nsr_partial2_theta2['2o-4n']),
        '2o-5n': coefficient * (compo_static_partial2_theta2 + compo_nsr_partial2_theta2['2o-5n'] +
                                compo_e_partial2_theta2['2o-5n'] + compo_e_nsr_partial2_theta2['2o-5n'])
        }

    potential_partial2_phi2 = {
        'n'    : coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['n'] +
                                compo_e_partial2_phi2['n'] + compo_e_nsr_partial2_phi2['n']),
        '2n'   : coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2n'] +
                                compo_e_partial2_phi2['2n'] + compo_e_nsr_partial2_phi2['2n']),
        '3n'   : coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['3n'] +
                                compo_e_partial2_phi2['3n'] + compo_e_nsr_partial2_phi2['3n']),
        '2o+n' : coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2o+n'] +
                                compo_e_partial2_phi2['2o+n'] + compo_e_nsr_partial2_phi2['2o+n']),
        '2o-n' : coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2o-n'] +
                                compo_e_partial2_phi2['2o-n'] + compo_e_nsr_partial2_phi2['2o-n']),
        '2o-2n': coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2o-2n'] +
                                compo_e_partial2_phi2['2o-2n'] + compo_e_nsr_partial2_phi2['2o-2n']),
        '2o-3n': coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2o-3n'] +
                                compo_e_partial2_phi2['2o-3n'] + compo_e_nsr_partial2_phi2['2o-3n']),
        '2o-4n': coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2o-4n'] +
                                compo_e_partial2_phi2['2o-4n'] + compo_e_nsr_partial2_phi2['2o-4n']),
        '2o-5n': coefficient * (compo_static_partial2_phi2 + compo_nsr_partial2_phi2['2o-5n'] +
                                compo_e_partial2_phi2['2o-5n'] + compo_e_nsr_partial2_phi2['2o-5n'])
        }

    potential_partial2_theta_phi = {
        'n'    : coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['n'] +
                                compo_e_partial2_theta_phi['n'] + compo_e_nsr_partial2_theta_phi['n']),
        '2n'   : coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2n'] +
                                compo_e_partial2_theta_phi['2n'] + compo_e_nsr_partial2_theta_phi['2n']),
        '3n'   : coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['3n'] +
                                compo_e_partial2_theta_phi['3n'] + compo_e_nsr_partial2_theta_phi['3n']),
        '2o+n' : coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2o+n'] +
                                compo_e_partial2_theta_phi['2o+n'] + compo_e_nsr_partial2_theta_phi['2o+n']),
        '2o-n' : coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2o-n'] +
                                compo_e_partial2_theta_phi['2o-n'] + compo_e_nsr_partial2_theta_phi['2o-n']),
        '2o-2n': coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2o-2n'] +
                                compo_e_partial2_theta_phi['2o-2n'] + compo_e_nsr_partial2_theta_phi['2o-2n']),
        '2o-3n': coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2o-3n'] +
                                compo_e_partial2_theta_phi['2o-3n'] + compo_e_nsr_partial2_theta_phi['2o-3n']),
        '2o-4n': coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2o-4n'] +
                                compo_e_partial2_theta_phi['2o-4n'] + compo_e_nsr_partial2_theta_phi['2o-4n']),
        '2o-5n': coefficient * (compo_static_partial2_theta_phi + compo_nsr_partial2_theta_phi['2o-5n'] +
                                compo_e_partial2_theta_phi['2o-5n'] + compo_e_nsr_partial2_theta_phi['2o-5n'])
        }

    return modes, potential, potential_partial_theta, potential_partial_phi, \
           potential_partial2_theta2, potential_partial2_phi2, potential_partial2_theta_phi
