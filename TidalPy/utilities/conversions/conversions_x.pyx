# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from TidalPy.exceptions import BadValueError

from libc.math cimport sqrt, cbrt

from TidalPy.constants cimport d_G, d_PI_DBL


cdef inline double cf_m2Au(double meters) noexcept nogil:

    return meters / 149597870700.0

cdef inline double cf_Au2m(double astronomical_units) noexcept nogil:

    return astronomical_units * 149597870700.0

cdef inline double cf_rads2days(double radians_per_second) noexcept nogil:

    return (2. * d_PI_DBL / radians_per_second) / 86400.

cdef inline double cf_days2rads(double days) noexcept nogil:

    return 2. * d_PI_DBL / (days * 86400.)

cdef inline double cf_sec2myr(double seconds) noexcept nogil:

    return seconds / 3.154e13

cdef inline double cf_myr2sec(double myrs) noexcept nogil:

    return myrs * 3.154e13

cdef inline double cf_orbital_motion2semi_a(
        double orbital_motion,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = d_G) noexcept nogil:

    cdef double semi_major_axis
    semi_major_axis = cbrt(
        G_to_use * (host_mass + target_mass) / (orbital_motion * orbital_motion)
        )

    return cbrt(G_to_use * (host_mass + target_mass) / (orbital_motion * orbital_motion))

cdef inline double cf_semi_a2orbital_motion(
        double semi_major_axis,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = d_G) noexcept nogil:

    cdef double orbital_motion
    orbital_motion = sqrt(G_to_use * (host_mass + target_mass) / (semi_major_axis * semi_major_axis * semi_major_axis))

    return orbital_motion

def m2Au(double meters):
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

    return cf_m2Au(meters)

def Au2m(double astronomical_units):
    """ Convert Meters to Astronomical Units

    Parameters
    ----------
    astronomical_units : double
        Distance in [Au]

    Returns
    -------
    meters : double
        Distance in [m]
    """

    return cf_Au2m(astronomical_units)

def rads2days(double radians_per_second):
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

    return cf_rads2days(radians_per_second)

def days2rads(double days):
    """ Convert from period [days] to frequency [rads s-1]

    Parameters
    ----------
    days : double
        Period in [days]

    Returns
    -------
    radians_per_second : double
        Frequency in [rads s-1]
    """

    return cf_days2rads(days)

def sec2myr(double seconds):
    """ Convert time from seconds to millions of years

    Parameters
    ----------
    seconds : double
        Time in [sec]

    Returns
    -------
    myrs : double
        Time in [Myr]
    """

    return cf_sec2myr(seconds)

def myr2sec(double myrs):
    """ Convert time from millions of years to seconds

    Parameters
    ----------
    myrs : double
        Time in [Myr]

    Returns
    -------
    seconds : double
        Time in [sec]
    """

    return cf_myr2sec(myrs)

def orbital_motion2semi_a(
        double orbital_motion,
        double host_mass,
        double target_mass = 0.0,
        double G_to_use = d_G):
    """ Convert orbital mean motion to semi-major axis (Kepler's 3rd law)

    Parameters
    ----------
    orbital_motion : double
        Orbital motion in [rads s-1]
    host_mass : double
        Central body's mass in [kg]
    target_mass : double, default = 0
        Target (or orbiting) body's mass in [kg]
    G_to_use : double, default = d_G
        Gravitational constant [N m2 kg-2]

    Returns
    -------
    semi_major_axis : double
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
        double G_to_use = d_G):
    """ Convert semi-major axis to mean orbital motion (Kepler's 3rd law)

    Parameters
    ----------
    semi_major_axis : double
        Semi-major axis in [m]
    host_mass : double
        Central body's mass in [kg]
    target_mass : double, default = 0
        Target (or orbiting) body's mass in [kg]
    G_to_use : double, default = d_G
        Gravitational constant [N m2 kg-2]

    Returns
    -------
    orbital_motion : double
        Orbital motion in [rads s-1]
    """

    if host_mass <= 0.:
        raise BadValueError('Host mass must be greater than zero.')
    if target_mass < 0.:
        raise BadValueError('Target mass must be greater than or equal to zero.')

    return cf_semi_a2orbital_motion(semi_major_axis, host_mass, target_mass, G_to_use)
