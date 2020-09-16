from typing import Tuple

import numpy as np

from .. import log
from ..tools.conversions import Au2m


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

    # Do a sanity check on the orbital motion / separation
    if orbital_period is not None:
        if orbital_freq is not None:
            log.info(f'Both orbital frequency and period were provided for {name}. Prioritizing frequency.')
            del orbital_period
        else:
            orbital_freq = np.asarray([2. * np.pi / (orbital_period * 24. * 60. * 60.)])

    if semi_major_axis is not None:
        semi_major_axis = np.asarray([semi_major_axis])
        if semi_major_axis_inau:
            semi_major_axis = Au2m(semi_major_axis)

    return orbital_freq, semi_major_axis, eccentricity
