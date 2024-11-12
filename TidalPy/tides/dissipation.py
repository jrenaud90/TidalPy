from typing import TYPE_CHECKING

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray


@njit(cacheable=True)
def calc_tidal_susceptibility(host_mass: float, target_radius: float, semi_major_axis: 'FloatArray') -> 'FloatArray':
    """ Calculate the tidal susceptibility for a given target radius, host mass, and their separation.

    Parameters
    ----------
    host_mass : float
        Mass of central host [kg]
    target_radius : float
        Radius of target body [m]
    semi_major_axis : FloatArray
        Semi-major axis [m]

    Returns
    -------
    tidal_susceptibility : FloatArray
        Tidal Susceptibility [N m]
    """

    tidal_susceptibility = (3. / 2.) * G * host_mass**2 * target_radius**5 / semi_major_axis**6

    return tidal_susceptibility


@njit(cacheable=True)
def calc_tidal_susceptibility_reduced(host_mass: float, target_radius: float) -> float:
    """ Calculate the tidal susceptibility (reduced) for a given target radius and host mass.

    The reduced tidal susceptibility excludes the semi-major axis.

    Parameters
    ----------
    host_mass : float
        Mass of central host [kg]
    target_radius : float
        Radius of target body [m]

    Returns
    -------
    tidal_susceptibility_reduced : np.ndarray
        Tidal Susceptibility [N m]
    """

    tidal_susceptibility_reduced = (3. / 2.) * G * host_mass**2 * target_radius**5

    return tidal_susceptibility_reduced
