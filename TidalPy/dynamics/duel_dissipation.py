import numpy as np
from scipy.constants import G

from ..performance import njit
from ..types import float_eps


# Duel Dissipation model's spin-rate derivative is the same as the single dissipation model

@njit
def semi_major_axis_derivative(semi_major_axis: np.ndarray, mass_1: float, mass_2: float,
                               spin_freq_1: np.ndarray, ztorque_1: np.ndarray, tidal_heating_1: np.ndarray,
                               spin_freq_2: np.ndarray, ztorque_2: np.ndarray, tidal_heating_2: np.ndarray
                               ) -> np.ndarray:
    """ Calculate the time derivative of the semi-major axis for a duel dissipating system

    See Joe Renaud's PhD Thesis (2019)

    Parameters
    ----------
    semi_major_axis : np.ndarray
        Semi-major axis in [m]
    mass_1 : float
        Mass of body 1 in [kg]
    mass_2 : float
        Mass of body 2 in [kg]
    spin_freq_1 : np.ndarray
        Spin frequency of body 1 in [rads s-1]
    ztorque_1 : np.ndarray
        Tidal polar torque of body 1 in [N m]
    tidal_heating_1 : np.ndarray
        Tidal heating of body 1 in [Watts]
    spin_freq_2 : np.ndarray
        Spin frequency of body 2 in [rads s-1]
    ztorque_2 : np.ndarray
        Tidal polar torque of body 2 in [N m]
    tidal_heating_2 : np.ndarray
        Tidal heating of body 2 in [Watts]

    Returns
    -------
    da_dt : np.ndarray
        Time derivative of the semi-major axis in [m s-1]
    """

    change_due_to_obj1 = spin_freq_1 * ztorque_1 + tidal_heating_1
    change_due_to_obj2 = spin_freq_2 * ztorque_2 + tidal_heating_2

    da_dt = (-2. * semi_major_axis**2 / (G * mass_1 * mass_2)) * (change_due_to_obj1 + change_due_to_obj2)

    return da_dt


@njit
def eccentricity_derivative(semi_major_axis: np.ndarray, eccentricity: np.ndarray, mass_1: float, mass_2: float,
                            spin_freq_1: np.ndarray, ztorque_1: np.ndarray, tidal_heating_1: np.ndarray,
                            spin_freq_2: np.ndarray, ztorque_2: np.ndarray, tidal_heating_2: np.ndarray) -> np.ndarray:
    """ Calculate the time derivative of the semi-major axis for a duel dissipating system

    See Joe Renaud's PhD Thesis (2019)

    Parameters
    ----------
    semi_major_axis : np.ndarray
        Semi-major axis in [m]
    eccentricity : np.ndarray
        Orbital Eccentricity
    mass_1 : float
        Mass of body 1 in [kg]
    mass_2 : float
        Mass of body 2 in [kg]
    spin_freq_1 : np.ndarray
        Spin frequency of body 1 in [rads s-1]
    ztorque_1 : np.ndarray
        Tidal polar torque of body 1 in [N m]
    tidal_heating_1 : np.ndarray
        Tidal heating of body 1 in [Watts]
    spin_freq_2 : np.ndarray
        Spin frequency of body 2 in [rads s-1]
    ztorque_2 : np.ndarray
        Tidal polar torque of body 2 in [N m]
    tidal_heating_2 : np.ndarray
        Tidal heating of body 2 in [Watts]

    Returns
    -------
    de_dt : np.ndarray
        Change in orbital eccentricity in [s-1]
    """

    # Check for bad values of eccentricity and set them to something workable for the equations
    bad_indices = eccentricity <= float_eps
    orbital_freq = np.sqrt(G * (mass_1 + mass_2) / semi_major_axis**3)
    e2_sqrt = np.sqrt(1 - eccentricity**2)

    change_due_to_obj1 = ztorque_1 * (1. - (spin_freq_1 / orbital_freq) * e2_sqrt) - \
                         (tidal_heating_1 / orbital_freq) * e2_sqrt
    change_due_to_obj2 = ztorque_2 * (1. - (spin_freq_2 / orbital_freq) * e2_sqrt) - \
                         (tidal_heating_2 / orbital_freq) * e2_sqrt

    de_dt = ((mass_1 + mass_2) * e2_sqrt / (mass_1 * mass_2 * semi_major_axis**2 * orbital_freq * eccentricity)) * \
            (change_due_to_obj1 + change_due_to_obj2)

    # The eccentricity is not going to change (unless perturbed) when e = 0 (what the bad_indices indicate)
    de_dt[bad_indices] = 0.

    return de_dt
