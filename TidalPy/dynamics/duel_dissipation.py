import numpy as np

from ..utilities.performance.numba import njit
from ..utilities.types import float_eps


# Duel Dissipation model's spin-rate derivative is the same as the single dissipation model

@njit
def semi_major_axis_derivative(semi_major_axis: np.ndarray, orbital_motion: np.ndarray,
                               mass_1: float, dU_dM_1: np.ndarray,
                               mass_2: float, dU_dM_2: np.ndarray) -> np.ndarray:
    """ Calculate the time derivative of the semi-major axis for a duel dissipating system
    See Boue and Efroimsky (2019, CMDA), Eq. 116
    Parameters
    ----------
    semi_major_axis : np.ndarray
        Semi-major axis in [m]
    orbital_motion : np.ndarray
        Orbital mean motion [rad s-1]
    mass_1 : float
        Mass of body 1 in [kg]
    dU_dM_1 : np.ndarray
        Derivative of body 1's tidal potential wrt mean anomaly
    mass_2 : float
        Mass of body 2 in [kg]
    dU_dM_2 : np.ndarray
        Derivative of body 2's tidal potential wrt mean anomaly
    Returns
    -------
    da_dt : np.ndarray
        Time derivative of the semi-major axis [m s-1]
    """

    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM_2 = -1. * beta_invr * mass_1 * dU_dM_2

    da_dt = (2. / (orbital_motion * semi_major_axis)) * (dR_dM_1 + dR_dM_2)

    return da_dt


@njit
def eccentricity_derivative(semi_major_axis: np.ndarray, orbital_motion: np.ndarray, eccentricity,
                            mass_1: float, dU_dM_1: np.ndarray, dU_dw_1: np.ndarray,
                            mass_2: float, dU_dM_2: np.ndarray, dU_dw_2: np.ndarray) -> np.ndarray:
    """ Calculate the time derivative of the eccentricity for a duel dissipating system
    See Boue and Efroimsky (2019, CMDA), Eq. 117
    Parameters
    ----------
    semi_major_axis : np.ndarray
        Semi-major axis in [m]
    orbital_motion : np.ndarray
        Orbital mean motion [rad s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    mass_1 : float
        Mass of body 1 in [kg]
    dU_dM_1 : np.ndarray
        Derivative of body 1's tidal potential wrt mean anomaly
    dU_dw_1 : np.ndarray
        Derivative of body 1's tidal potential wrt pericentre
    mass_2 : float
        Mass of body 2 in [kg]
    dU_dM_2 : np.ndarray
        Derivative of body 2's tidal potential wrt mean anomaly
    dU_dw_2 : np.ndarray
        Derivative of body 2's tidal potential wrt pericentre
    Returns
    -------
    da_dt : np.ndarray
        Time derivative of the eccentricity [s-1]
    """
    e_term1 = np.sqrt(1. - eccentricity * eccentricity)
    denom = orbital_motion * semi_major_axis * semi_major_axis * eccentricity


    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM_2 = -1. * beta_invr * mass_1 * dU_dM_2
    dR_dM = dR_dM_1 + dR_dM_2

    dR_dw_1 = -1. * beta_invr * mass_2 * dU_dw_1
    dR_dw_2 = -1. * beta_invr * mass_1 * dU_dw_2
    de_dt = (e_term1 / denom) * (e_term1 * dR_dM - (dR_dw_1 + dR_dw_2))

    # Correct for zero eccentricity
    de_dt[np.abs(denom) <= float_eps] = 0.

    return de_dt