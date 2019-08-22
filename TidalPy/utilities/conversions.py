from typing import Tuple

import numpy as np

from ..constants import G
from ..performance import njit
from ..types import FloatArray


@njit
def m2Au(meters: FloatArray) -> FloatArray:
    """ Convert Meters to Astronomical Units

    Parameters
    ----------
    meters : FloatArray
        Distance in [m]

    Returns
    -------
    astronomical_units : FloatArray
        Distance in [Au]
    """

    astronomical_units = meters / 1.496e11

    return astronomical_units


@njit
def Au2m(astronomical_units: FloatArray) -> FloatArray:
    """ Convert Meters to Astronomical Units

    Parameters
    ----------
    astronomical_units : FloatArray
        Distance in [Au]

    Returns
    -------
    meters : FloatArray
        Distance in [m]
    """

    meters = astronomical_units * 1.496e11

    return meters


@njit
def rads2days(radians_per_second: FloatArray) -> FloatArray:
    """ Convert from frequency [rads s-1] to period [days]

    Parameters
    ----------
    radians_per_second : FloatArray
        Frequency in [rads s-1]

    Returns
    -------
    days : FloatArray
        Period in [days]
    """

    days = (2. * np.pi / radians_per_second) / 86400.

    return days


@njit
def days2rads(days: FloatArray) -> FloatArray:
    """ Convert from period [days] to frequency [rads s-1]

    Parameters
    ----------
    days : FloatArray
        Period in [days]

    Returns
    -------
    radians_per_second : FloatArray
        Frequency in [rads s-1]
    """

    radians_per_second = (2. * np.pi / (days * 86400.))

    return radians_per_second


@njit
def orbital_motion2semi_a(orbital_motion: FloatArray, host_mass: float, target_mass: float = 0.) -> FloatArray:
    """ Convert orbital mean motion to semi-major axis (Kepler's 3rd law)

    Parameters
    ----------
    orbital_motion : FloatArray
        Orbital motion in [rads s-1]
    host_mass : float
        Central body's mass in [kg]
    target_mass : float = 0.
        Target (or orbiting) body's mass in [kg]

    Returns
    -------
    semi_major_axis : FloatArray
        Semi-major axis in [m]
    """

    semi_major_axis = (G * (host_mass + target_mass) / orbital_motion**2)**(1 / 3)

    return semi_major_axis


@njit
def semi_a2orbital_motion(semi_major_axis: FloatArray, host_mass: float, target_mass: float = 0.) -> FloatArray:
    """ Convert semi-major axis to mean orbital motion (Kepler's 3rd law)

    Parameters
    ----------
    semi_major_axis : FloatArray
        Semi-major axis in [m]
    host_mass : float
        Central body's mass in [kg]
    target_mass : float = 0.
        Target (or orbiting) body's mass in [kg]

    Returns
    -------
    orbital_motion : FloatArray
        Orbital motion in [rads s-1]
    """

    orbital_motion = (G * (host_mass + target_mass) / semi_major_axis**3)**(1 / 2)

    return orbital_motion


@njit
def sec2myr(seconds: FloatArray) -> FloatArray:
    """ Convert time from seconds to millions of years

    Parameters
    ----------
    seconds : FloatArray
        Time in [sec]

    Returns
    -------
    myrs : FloatArray
        Time in [Myr]
    """

    return seconds / 3.154e13


@njit
def myr2sec(myrs: FloatArray) -> FloatArray:
    """ Convert time from millions of years to seconds

    Parameters
    ----------
    myrs : FloatArray
        Time in [Myr]

    Returns
    -------
    seconds : FloatArray
        Time in [sec]
    """

    return myrs * 3.154e13


@njit
def convert_to_hms(seconds: float) -> Tuple[int, int, int, float]:
    """ Convert seconds to a tuple of days, hours, minutes, seconds

    Parameters
    ----------
    seconds : float
        Time in seconds

    Returns
    -------
    days : int
        Days
    hours : int
        Hours
    minutes : int
        Minutes
    seconds : float
        Remaining seconds after conversion
    """

    days = int(seconds / (24. * 3600.))
    seconds = seconds % (24. * 3600.)
    hours = int(seconds / 3600.)
    seconds = seconds % 3600.
    minutes = int(seconds / 60.)
    seconds = seconds % 60.

    return days, hours, minutes, seconds
