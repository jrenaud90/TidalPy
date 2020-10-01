from typing import Tuple

import numpy as np

from ..utilities.performance.numba import njit
from ..utilities.types import float_eps, FloatArray


# Duel Dissipation model's spin-rate derivative is the same as the single dissipation model

@njit(cacheable=True)
def semi_major_axis_derivative(semi_major_axis: FloatArray, orbital_motion: FloatArray,
                               mass_1: float, dU_dM_1: FloatArray,
                               mass_2: float, dU_dM_2: FloatArray) -> FloatArray:
    """ Calculate the time derivative of the semi-major axis for a dual dissipating system

    See Boue and Efroimsky (2019, CMDA), Eq. 116

    Parameters
    ----------
    semi_major_axis : FloatArray
        Semi-major axis in [m]
    orbital_motion : FloatArray
        Orbital mean motion [rad s-1]
    mass_1 : float
        Mass of body 1 in [kg]
    dU_dM_1 : FloatArray
        Derivative of body 1's tidal potential wrt mean anomaly
    mass_2 : float
        Mass of body 2 in [kg]
    dU_dM_2 : FloatArray
        Derivative of body 2's tidal potential wrt mean anomaly

    Returns
    -------
    da_dt : FloatArray
        Time derivative of the semi-major axis [m s-1]
    """

    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM_2 = -1. * beta_invr * mass_1 * dU_dM_2
    dR_dM = dR_dM_1 + dR_dM_2

    da_dt = (2. / (orbital_motion * semi_major_axis)) * dR_dM

    return da_dt

@njit(cacheable=True)
def eccentricity_derivative(semi_major_axis: float, orbital_motion: float, eccentricity: float,
                            mass_1: float, dU_dM_1: float, dU_dw_1: float,
                            mass_2: float, dU_dM_2: float, dU_dw_2: float) -> float:
    """ Calculate the time derivative of the eccentricity for a dual dissipating system - NonArrays Only

    See Boue and Efroimsky (2019, CMDA), Eq. 117

    Parameters
    ----------
    semi_major_axis : float
        Semi-major axis in [m]
    orbital_motion : float
        Orbital mean motion [rad s-1]
    eccentricity : float
        Orbital eccentricity
    mass_1 : float
        Mass of body 1 in [kg]
    dU_dM_1 : float
        Derivative of body 1's tidal potential wrt mean anomaly
    dU_dw_1 : float
        Derivative of body 1's tidal potential wrt pericentre
    mass_2 : float
        Mass of body 2 in [kg]
    dU_dM_2 : float
        Derivative of body 2's tidal potential wrt mean anomaly
    dU_dw_2 : float
        Derivative of body 2's tidal potential wrt pericentre

    Returns
    -------
    de_dt : float
        Time derivative of the eccentricity [s-1]
    """

    e_term1 = np.sqrt(1. - eccentricity * eccentricity)
    denom = orbital_motion * semi_major_axis * semi_major_axis * eccentricity

    if abs(denom) <= float_eps:
        # Correct for zero eccentricity
        return 0.

    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM_2 = -1. * beta_invr * mass_1 * dU_dM_2
    dR_dM = dR_dM_1 + dR_dM_2

    dR_dw_1 = -1. * beta_invr * mass_2 * dU_dw_1
    dR_dw_2 = -1. * beta_invr * mass_1 * dU_dw_2
    de_dt = (e_term1 / denom) * (e_term1 * dR_dM - (dR_dw_1 + dR_dw_2))

    return de_dt


@njit(cacheable=True)
def eccentricity_derivative_array(semi_major_axis: FloatArray, orbital_motion: FloatArray, eccentricity,
                                  mass_1: float, dU_dM_1: FloatArray, dU_dw_1: FloatArray,
                                  mass_2: float, dU_dM_2: FloatArray, dU_dw_2: FloatArray) -> FloatArray:
    """ Calculate the time derivative of the eccentricity for a dual dissipating system - Arrays Only

    See Boue and Efroimsky (2019, CMDA), Eq. 117

    Parameters
    ----------
    semi_major_axis : FloatArray
        Semi-major axis in [m]
    orbital_motion : FloatArray
        Orbital mean motion [rad s-1]
    eccentricity : FloatArray
        Orbital eccentricity
    mass_1 : float
        Mass of body 1 in [kg]
    dU_dM_1 : FloatArray
        Derivative of body 1's tidal potential wrt mean anomaly
    dU_dw_1 : FloatArray
        Derivative of body 1's tidal potential wrt pericentre
    mass_2 : float
        Mass of body 2 in [kg]
    dU_dM_2 : FloatArray
        Derivative of body 2's tidal potential wrt mean anomaly
    dU_dw_2 : FloatArray
        Derivative of body 2's tidal potential wrt pericentre

    Returns
    -------
    de_dt : FloatArray
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
    de_dt[np.abs(denom * dR_dM) <= float_eps] = 0.

    return de_dt

@njit(cacheable=True)
def semia_eccen_derivatives(semi_major_axis: float, orbital_motion: float, eccentricity: float,
                            mass_1: float, dU_dM_1: float, dU_dw_1: float,
                            mass_2: float, dU_dM_2: float, dU_dw_2: float) \
        -> Tuple[float, float]:
    """ Calculate the time derivatives of semi-major axis and eccentricity for a dual dissipating system - NonArrays Only

    See Boue and Efroimsky (2019, CMDA), Eqs. 116 and 117

    Parameters
    ----------
    semi_major_axis : float
        Semi-major axis in [m]
    orbital_motion : float
        Orbital mean motion [rad s-1]
    eccentricity : float
        Orbital eccentricity
    mass_1 : float
        Mass of body 1 in [kg]
    dU_dM_1 : float
        Derivative of body 1's tidal potential wrt mean anomaly
    dU_dw_1 : float
        Derivative of body 1's tidal potential wrt pericentre
    mass_2 : float
        Mass of body 2 in [kg]
    dU_dM_2 : float
        Derivative of body 2's tidal potential wrt mean anomaly
    dU_dw_2 : float
        Derivative of body 2's tidal potential wrt pericentre

    Returns
    -------
    da_dt : float
        Time derivative of the semi-major axis [m s-1]
    de_dt : float
        Time derivative of the eccentricity [s-1]
    """

    # Calculate semi-major axis derivative
    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM_2 = -1. * beta_invr * mass_1 * dU_dM_2
    dR_dM = dR_dM_1 + dR_dM_2

    da_dt = (2. / (orbital_motion * semi_major_axis)) * dR_dM

    # Calculate eccentricity derivative
    e_term1 = np.sqrt(1. - eccentricity * eccentricity)
    denom = orbital_motion * semi_major_axis * semi_major_axis * eccentricity

    if abs(denom) <= float_eps:
        # Correct for zero eccentricity
        de_dt = 0.
    else:
        dR_dw_1 = -1. * beta_invr * mass_2 * dU_dw_1
        dR_dw_2 = -1. * beta_invr * mass_1 * dU_dw_2
        de_dt = (e_term1 / denom) * (e_term1 * dR_dM - (dR_dw_1 + dR_dw_2))

    return da_dt, de_dt

@njit(cacheable=True)
def semia_eccen_derivatives_array(semi_major_axis: np.ndarray, orbital_motion: np.ndarray, eccentricity: np.ndarray,
                                  mass_1: float, dU_dM_1: np.ndarray, dU_dw_1: np.ndarray,
                                  mass_2: float, dU_dM_2: np.ndarray, dU_dw_2: np.ndarray) \
        -> Tuple[np.ndarray, np.ndarray]:
    """ Calculate the time derivatives of semi-major axis and eccentricity for a dual dissipating system - Arrays Only

    See Boue and Efroimsky (2019, CMDA), Eqs. 116 and 117

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
    da_dt : FloatArray
        Time derivative of the semi-major axis [m s-1]
    de_dt : FloatArray
        Time derivative of the eccentricity [s-1]
    """

    # Calculate semi-major axis derivative
    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM_2 = -1. * beta_invr * mass_1 * dU_dM_2
    dR_dM = dR_dM_1 + dR_dM_2

    da_dt = (2. / (orbital_motion * semi_major_axis)) * dR_dM

    # Calculate eccentricity derivative
    e_term1 = np.sqrt(1. - eccentricity * eccentricity)
    denom = orbital_motion * semi_major_axis * semi_major_axis * eccentricity

    dR_dw_1 = -1. * beta_invr * mass_2 * dU_dw_1
    dR_dw_2 = -1. * beta_invr * mass_1 * dU_dw_2
    de_dt = (e_term1 / denom) * (e_term1 * dR_dM - (dR_dw_1 + dR_dw_2))

    # Correct for zero eccentricity
    de_dt[np.abs(denom * dR_dM) <= float_eps] = 0.

    return da_dt, de_dt

