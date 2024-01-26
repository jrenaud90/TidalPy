# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from TidalPy.exceptions import BadValueError

from libc.math cimport sqrt, cbrt, M_PI

from TidalPy.utilities.constants cimport G


@njit(cacheable=True)
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


@njit(cacheable=True)
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


@njit(cacheable=True)
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


@njit(cacheable=True)
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


@njit(cacheable=True)
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


@njit(cacheable=True)
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

cdef double cf_orbital_motion2semi_a(
        double orbital_motion,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = G):

    cdef double semi_major_axis
    semi_major_axis = cbrt(
        G_to_use * (host_mass + target_mass) / (orbital_motion * orbital_motion)
        )

    return semi_major_axis


cdef double cf_semi_a2orbital_motion(
        double semi_major_axis,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = G):

    cdef double orbital_motion
    orbital_motion = sqrt(
        G_to_use * (host_mass + target_mass) / (semi_major_axis * semi_major_axis * semi_major_axis)
        )

    return orbital_motion

def orbital_motion2semi_a(
        double orbital_motion,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = G):
    """ Convert orbital mean motion to semi-major axis (Kepler's 3rd law)

    Parameters
    ----------
    orbital_motion : FloatArray
        Orbital motion in [rads s-1]
    host_mass : float
        Central body's mass in [kg]
    target_mass : double, default = 0
        Target (or orbiting) body's mass in [kg]
    G_to_use : double, default = G
        Gravitational constant [N m2 kg-2]

    Returns
    -------
    semi_major_axis : FloatArray
        Semi-major axis in [m]
    """

    if host_mass <= 0.:
        raise BadValueError('Host mass must be greater than zero.')
    if target_mass < 0.:
        raise BadValueError('Target mass must be greater than or equal to zero.')

    return cf_orbital_motion2semi_a(orbital_motion, host_mass, target_mass, G_to_use)

def semi_a2orbital_motion(
        double semi_major_axis,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = G):
    """ Convert semi-major axis to mean orbital motion (Kepler's 3rd law)

    Parameters
    ----------
    semi_major_axis : double
        Semi-major axis in [m]
    host_mass : double
        Central body's mass in [kg]
    target_mass : double, default = 0
        Target (or orbiting) body's mass in [kg]
    G_to_use : double, default = G
        Gravitational constant [N m2 kg-2]

    Returns
    -------
    orbital_motion : FloatArray
        Orbital motion in [rads s-1]
    """

    if host_mass <= 0.:
        raise BadValueError('Host mass must be greater than zero.')
    if target_mass < 0.:
        raise BadValueError('Target mass must be greater than or equal to zero.')

    return cf_semi_a2orbital_motion(semi_major_axis, host_mass, target_mass, G_to_use)
