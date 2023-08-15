""" Functions to help calculate the mass and gravity of a spherical body
"""
from typing import Tuple

import numpy as np

from TidalPy.constants import G

from TidalPy.utilities.performance import njit


@njit(cacheable=True)
def calculate_mass_gravity_arrays(radius_array: np.ndarray, density_array: np.ndarray,
                                  gravity_constant: float = G) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Calculate the volume, mass, and gravity arrays from a radius and density array.

    Parameters
    ----------
    radius_array : np.ndarray
        Radius array throughout planet [m or units that match `gravity_constant`]
        This does not need to be evenly spaced.
    density_array : np.ndarray
        Density array throughout the planet [kg m-3 or units that match `gravity_constant`]
        Each element of the array must match the density at the respective elements in `radius_array`
    gravity_constant : float = G (from TidalPy.constants)
        The user can provide an alternative G if they are using a different unit system.

    Returns
    -------
    output : Tuple[np.ndarray, np.ndarray, np.ndarray]
        volume_array : np.ndarray
            Volume of each spherical shell [m3]
            Defined at each element of the `radius_array`
        mass_array : np.ndarray
            Mass of each spherical shell [kg]
            Defined at each element of the `radius_array`
        gravity_array : np.ndarray
            Acceleration due to gravity at the top of each spherical shell [m s-3]
            Defined at each element of the `radius_array`

    """

    num_shells = len(radius_array)

    # Optimizations
    r2 = radius_array * radius_array
    r3 = r2 * radius_array

    # Calculate volume and mass
    volume_base_shell = np.asarray(((4. / 3.) * np.pi * r3[0],))
    volume_other_shells = (4. / 3.) * np.pi * (r3[1:] - r3[:-1])
    volume_array = np.concatenate((volume_base_shell, volume_other_shells))
    mass_array = volume_array * density_array

    # Calculate the mass below each spherical shell. Assume that the mass below a shell also includes the mass of that
    #    shell.
    mass_below_array = np.asarray([np.sum(mass_array[0:i+1]) for i in range(num_shells)])

    # Calculate the acceleration due to gravity at the surface of each shell
    gravity_array = gravity_constant * mass_below_array / r2

    return volume_array, mass_array, gravity_array
