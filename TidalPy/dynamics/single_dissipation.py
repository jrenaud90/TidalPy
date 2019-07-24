import numpy as np
from scipy.constants import G

from ..performance import njit
from ..types import float_eps


@njit
def spin_rate_derivative(ztorque: np.ndarray, moment_of_inertia: float) -> np.ndarray:
    """ Calculate the time derivative of the spin frequency for a single-body or duel dissipation system

    See Ferraz-Mello et. al. (2008)

    Parameters
    ----------
    ztorque : np.ndarray
        Tidal polar torque in [N m]
    moment_of_inertia : float
        Planet's moment of inertia in [kg m2]

    Returns
    -------
    dspin_dt : np.ndarray
        Change of spin frequency in [rad s-2]
    """

    dspin_dt = ztorque / moment_of_inertia

    return dspin_dt


@njit
def semi_major_axis_derivative(semi_major_axis: np.ndarray, mass_host: float, mass_target: float,
                               spin_freq: np.ndarray, ztorque: np.ndarray, tidal_heating: np.ndarray) -> np.ndarray:
    """ Calculate the time derivative of the semi-major axis for a single-body system

    See Joe Renaud's PhD Thesis (2019)

    Parameters
    ----------
    semi_major_axis : np.ndarray
        Semi-major axis in [m]
    mass_host : float
        Mass of tidal host in [kg]
    mass_target : float
        Mass of target body in [kg]
    spin_freq : np.ndarray
        Spin frequency of target body in [rads s-1]
    ztorque : np.ndarray
        Tidal polar torque of target body in [N m]
    tidal_heating : np.ndarray
        Tidal heating of target body in [Watts]


    Returns
    -------
    da_dt : np.ndarray
        Time derivative of the semi-major axis in [m s-1]
    """

    change_due_to_obj1 = spin_freq * ztorque + tidal_heating
    change_due_to_obj2 = spin_freq * ztorque + tidal_heating

    da_dt = (-2. * semi_major_axis**2 / (G * mass_host * mass_target)) * (change_due_to_obj1 + change_due_to_obj2)

    return da_dt


@njit
def eccentricity_derivative(semi_major_axis: np.ndarray, eccentricity: np.ndarray, mass_host: float, mass_target: float,
                            spin_freq: np.ndarray, ztorque: np.ndarray, tidal_heating: np.ndarray) -> np.ndarray:
    """ Calculate the time derivative of the semi-major axis for a single-body system

    See Joe Renaud's PhD Thesis (2019)

    Parameters
    ----------
    semi_major_axis : np.ndarray
        Semi-major axis in [m]
    eccentricity : np.ndarray
        Orbital Eccentricity
    mass_host : float
        Mass of tidal host in [kg]
    mass_target : float
        Mass of target body in [kg]
    spin_freq : np.ndarray
        Spin frequency of target body in [rads s-1]
    ztorque : np.ndarray
        Tidal polar torque of target body in [N m]
    tidal_heating : np.ndarray
        Tidal heating of target body in [Watts]

    Returns
    -------
    de_dt : np.ndarray
        Change in orbital eccentricity in [s-1]
    """

    # Check for bad values of eccentricity and set them to something workable for the equations
    bad_indices = eccentricity <= float_eps
    orbital_freq = np.sqrt(G * (mass_host + mass_target) / semi_major_axis**3)
    e2_sqrt = np.sqrt(1 - eccentricity**2)

    change_due_to_target = ztorque * (1. - (spin_freq / orbital_freq) * e2_sqrt) - \
                           (tidal_heating / orbital_freq) * e2_sqrt

    de_dt = ((mass_host + mass_target) * e2_sqrt /
             (mass_host * mass_target * semi_major_axis**2 * orbital_freq * eccentricity)) * change_due_to_target

    # The eccentricity is not going to change (unless perturbed) when e = 0 (what the bad_indices indicate)
    de_dt[bad_indices] = 0.

    return de_dt
