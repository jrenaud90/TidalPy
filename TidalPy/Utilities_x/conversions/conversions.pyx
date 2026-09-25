# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Common unit and orbital-element conversions.

The Kepler conversions are the C++ ``c_semi_a2orbital_motion`` and ``c_orbital_motion2semi_a`` (``conversions_.hpp``),
which C++ callers use directly; these functions check the inputs and form the gravitational parameter.
"""

from TidalPy.constants cimport (
    d_PI, d_SECONDS_PER_MYR, tidalpy_config_ptr, get_shared_config_address, set_tidalpy_config_ptr)

# Wire this DLL's shared pointer to the process-wide TidalPy config singleton.
set_tidalpy_config_ptr(get_shared_config_address())


def m2Au(double meters):
    """ Convert Meters to Astronomical Units

    Parameters
    ----------
    meters : float
        Distance in [m]

    Returns
    -------
    astronomical_units : float
        Distance in [Au]
    """

    return meters / tidalpy_config_ptr.d_AU

def Au2m(double astronomical_units):
    """ Convert Astronomical Units to Meters

    Parameters
    ----------
    astronomical_units : float
        Distance in [Au]

    Returns
    -------
    meters : float
        Distance in [m]
    """

    return astronomical_units * tidalpy_config_ptr.d_AU

def rads2days(double radians_per_second):
    """ Convert from frequency [rads s-1] to period [days]

    Parameters
    ----------
    radians_per_second : float
        Frequency in [rads s-1]

    Returns
    -------
    days : float
        Period in [days]
    """

    return (2. * d_PI / radians_per_second) / 86400.

def days2rads(double days):
    """ Convert from period [days] to frequency [rads s-1]

    Parameters
    ----------
    days : float
        Period in [days]

    Returns
    -------
    radians_per_second : float
        Frequency in [rads s-1]
    """

    return 2. * d_PI / (days * 86400.)

def sec2myr(double seconds):
    """ Convert time from seconds to millions of Julian years (365.25 days each)

    Parameters
    ----------
    seconds : float
        Time in [sec]

    Returns
    -------
    myrs : float
        Time in [Myr]
    """

    return seconds / d_SECONDS_PER_MYR

def myr2sec(double myrs):
    """ Convert time from millions of Julian years (365.25 days each) to seconds

    Parameters
    ----------
    myrs : float
        Time in [Myr]

    Returns
    -------
    seconds : float
        Time in [sec]
    """

    return myrs * d_SECONDS_PER_MYR

def orbital_motion2semi_a(
        double orbital_motion,
        double host_mass,
        double target_mass = 0.0,
        G_to_use = None):
    """ Convert orbital mean motion to semi-major axis (Kepler's 3rd law)

    Parameters
    ----------
    orbital_motion : float
        Orbital motion in [rads s-1]
    host_mass : float
        Central body's mass in [kg]
    target_mass : float, default = 0
        Target (or orbiting) body's mass in [kg]
    G_to_use : float, optional
        Gravitational constant [m3 kg-1 s-2]. ``None`` (the default) reads the TidalPy config's value at call
        time.

    Returns
    -------
    semi_major_axis : float
        Semi-major axis in [m]

    Raises
    ------
    ValueError
        If the orbital motion or the host mass is not positive, the target mass is negative, or ``G_to_use``
        is not positive.

    Assumptions
    -----------
    Two-body Keplerian orbit: n^2 a^3 = G (M_host + M_target).
    """

    # Read at call time so a reinitialized config is honored.
    cdef double G_value = tidalpy_config_ptr.d_G if G_to_use is None else <double>G_to_use

    # Written as `not (x > 0)` so a NaN is refused as well.
    if not (orbital_motion > 0.):
        raise ValueError(f'Orbital motion must be greater than zero; got {orbital_motion} rad s-1.')
    if host_mass <= 0.:
        raise ValueError('Host mass must be greater than zero.')
    if target_mass < 0.:
        raise ValueError('Target mass must be greater than or equal to zero.')
    if not (G_value > 0.):
        raise ValueError(f'G_to_use must be greater than zero; got {G_value}.')

    return c_orbital_motion2semi_a(orbital_motion, G_value * (host_mass + target_mass))

def semi_a2orbital_motion(
        double semi_major_axis,
        double host_mass,
        double target_mass = 0.0,
        G_to_use = None):
    """ Convert semi-major axis to mean orbital motion (Kepler's 3rd law)

    Parameters
    ----------
    semi_major_axis : float
        Semi-major axis in [m]
    host_mass : float
        Central body's mass in [kg]
    target_mass : float, default = 0
        Target (or orbiting) body's mass in [kg]
    G_to_use : float, optional
        Gravitational constant [m3 kg-1 s-2]. ``None`` (the default) reads the TidalPy config's value at call
        time.

    Returns
    -------
    orbital_motion : float
        Orbital motion in [rads s-1]

    Raises
    ------
    ValueError
        If the semi-major axis or the host mass is not positive, the target mass is negative, or ``G_to_use``
        is not positive.

    Assumptions
    -----------
    Two-body Keplerian orbit: n^2 a^3 = G (M_host + M_target).
    """

    # Read at call time so a reinitialized config is honored.
    cdef double G_value = tidalpy_config_ptr.d_G if G_to_use is None else <double>G_to_use

    # Written as `not (x > 0)` so a NaN is refused as well.
    if not (semi_major_axis > 0.):
        raise ValueError(f'Semi-major axis must be greater than zero; got {semi_major_axis} m.')
    if host_mass <= 0.:
        raise ValueError('Host mass must be greater than zero.')
    if target_mass < 0.:
        raise ValueError('Target mass must be greater than or equal to zero.')
    if not (G_value > 0.):
        raise ValueError(f'G_to_use must be greater than zero; got {G_value}.')

    return c_semi_a2orbital_motion(semi_major_axis, G_value * (host_mass + target_mass))
