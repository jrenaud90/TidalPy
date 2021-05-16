from typing import Tuple

import numpy as np

from ... import log
from ...toolbox.conversions import Au2m, days2rads, rads2days


def pull_out_orbit_from_config(world_config: dict) -> Tuple[float, float, float]:
    """ Parse a user-provided configuration dict and pull out any orbital information.

    Also perform some sanity checks and conversions.

    Parameters
    ----------
    world_config : dict
        User-provided world configuration

    Returns
    -------
    orbital_freq : float
        Orbital mean motion [rad s-1]
    semi_major_axis : float
        Orbital separation [m]
    eccentricity : float
        Orbital eccentricity
    """

    # See what is in the configuration file, if anything.
    orbital_freq = None
    if 'orbital_freq' in world_config:
        orbital_freq = world_config['orbital_freq']
    elif 'orbital_frequency' in world_config:
        orbital_freq = world_config['orbital_frequency']
    elif 'orbital_motion' in world_config:
        orbital_freq = world_config['orbital_motion']
    elif 'orbital_mean_motion' in world_config:
        orbital_freq = world_config['orbital_mean_motion']

    orbital_period = world_config.get('orbital_period', None)
    semi_major_axis = world_config.get('semi_major_axis', None)
    semi_major_axis_inau = world_config.get('semi_major_axis_in_au', False)
    eccentricity = world_config.get('eccentricity', None)
    name = world_config['name']

    if orbital_freq is not None:
        if np.any(orbital_freq > 0.001):
            log.warning(f'Unusual Value Encountered in {name}: User-provided orbital frequency set in [rad s-1] but '
                        f'value is very large, {np.max(orbital_freq):0.3E} '
                        f'(equivalent to {rads2days(np.max(orbital_freq)):0.2f} day period).')
        elif np.any(orbital_freq < 1.e-7):
            log.warning(f'Unusual Value Encountered in {name}: User-provided orbital frequency set in [rad s-1] but '
                        f'value is very small, {np.max(orbital_freq):0.3E} '
                        f'(equivalent to {rads2days(np.max(orbital_freq)):0.2f} day period).')

    # Do a sanity check on the orbital motion / separation
    if orbital_period is not None:
        if orbital_freq is not None:
            log.warning(f'Both orbital frequency and period were provided for {name}. Prioritizing frequency.')
            del orbital_period
        else:
            if np.any(orbital_period > 500.):
                log.warning(f'Unusual Value Encountered in {name}: User-provided orbital period set in [days] but '
                            f'value is very large, {np.max(orbital_period):0.2f}.')
            elif np.any(orbital_period < .05):
                log.warning(f'Unusual Value Encountered in {name}: User-provided orbital period set in [days] but '
                            f'value is very small, {np.min(orbital_period):0.4f}.')
            orbital_freq = days2rads(orbital_period)

    if semi_major_axis is not None:
        semi_major_axis = semi_major_axis
        if semi_major_axis_inau:
            if np.any(semi_major_axis > 1000.):
                log.warning(f'Unusual Value Encountered in {name}: User-provided semi-major axis set in [Au] but value '
                            f'is very large, {np.max(semi_major_axis)}.')
            semi_major_axis = Au2m(semi_major_axis)
        else:
            if np.any(semi_major_axis < 1.e6):
                log.warning(f'Unusual Value Encountered in {name}: User-provided semi-major axis set in [m] but value '
                            f'is very small, {np.min(semi_major_axis)}.')



    return orbital_freq, semi_major_axis, eccentricity
