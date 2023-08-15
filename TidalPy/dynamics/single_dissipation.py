from typing import Tuple, TYPE_CHECKING

import numpy as np

from TidalPy.utilities.types import float_eps
from TidalPy.utilities.performance.numba import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray


@njit(cacheable=True)
def spin_rate_derivative(dU_dO: 'FloatArray', moment_of_inertia: float, host_mass: float) -> 'FloatArray':
    """ Calculate the time derivative of the spin frequency for a single-body or single or dual dissipation system

    See Ferraz-Mello et. al. (2008)

    Parameters
    ----------
    dU_dO : FloatArray
        Partial derivative of the tidal potential wrt the target planet's orbital node
    moment_of_inertia : float
        Planet's moment of inertia in [kg m2]
    host_mass : float
        Mass of the tidal host [kg]

    Returns
    -------
    dspin_dt : np.ndarray
        Change of spin frequency in [rad s-2]
    """

    dspin_dt = (host_mass / moment_of_inertia) * dU_dO

    return dspin_dt


@njit(cacheable=True)
def semi_major_axis_derivative(
    semi_major_axis: 'FloatArray', orbital_motion: 'FloatArray',
    mass_1: float, dU_dM_1: 'FloatArray', mass_2: float
    ) -> 'FloatArray':
    """ Calculate the time derivative of the semi-major axis for a signle-body dissipating system

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

    Returns
    -------
    da_dt : FloatArray
        Time derivative of the semi-major axis [m s-1]
    """

    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM = dR_dM_1

    da_dt = (2. / (orbital_motion * semi_major_axis)) * dR_dM

    return da_dt


@njit(cacheable=True)
def eccentricity_derivative(
    semi_major_axis: 'FloatArray', orbital_motion: 'FloatArray', eccentricity: 'FloatArray',
    mass_1: float, dU_dM_1: 'FloatArray', dU_dw_1: 'FloatArray', mass_2: float
    ) -> 'FloatArray':
    """ Calculate the time derivative of the eccentricity for a single-body dissipating system

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

    Returns
    -------
    de_dt : FloatArray
        Time derivative of the eccentricity [s-1]
    """

    e_term1 = np.sqrt(1. - eccentricity * eccentricity)
    denom = orbital_motion * semi_major_axis * semi_major_axis * eccentricity

    beta_invr = (mass_1 + mass_2) / (mass_1 * mass_2)
    dR_dM_1 = -1. * beta_invr * mass_2 * dU_dM_1
    dR_dM = dR_dM_1

    dR_dw_1 = -1. * beta_invr * mass_2 * dU_dw_1

    # Correct for zero eccentricity
    de_dt = (np.abs(denom) <= float_eps) * 0. + \
            (np.abs(denom) > float_eps) * (e_term1 / denom) * (e_term1 * dR_dM - dR_dw_1)

    return de_dt


@njit(cacheable=True)
def semia_eccen_derivatives(
    semi_major_axis: 'FloatArray', orbital_motion: 'FloatArray', eccentricity: 'FloatArray',
    mass_1: float, dU_dM_1: 'FloatArray', dU_dw_1: 'FloatArray', mass_2: float
    ) -> Tuple['FloatArray', 'FloatArray']:
    """ Calculate the time derivatives of semi-major axis and eccentricity for a single-body dissipating system

    See Boue and Efroimsky (2019, CMDA), Eqs. 116 and 117

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
    dR_dM = dR_dM_1

    dR_dw_1 = -1. * beta_invr * mass_2 * dU_dw_1

    da_dt = (2. / (orbital_motion * semi_major_axis)) * dR_dM

    # Calculate eccentricity derivative
    e_term1 = np.sqrt(1. - eccentricity * eccentricity)
    denom = orbital_motion * semi_major_axis * semi_major_axis * eccentricity

    # Correct for zero eccentricity
    de_dt = (np.abs(denom) <= float_eps) * 0. + \
            (np.abs(denom) > float_eps) * (e_term1 / denom) * (e_term1 * dR_dM - dR_dw_1)

    return da_dt, de_dt
