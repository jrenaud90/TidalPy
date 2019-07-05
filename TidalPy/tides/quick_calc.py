from typing import Tuple, List, Union

import numpy as np

from ..performance import njit
from ..types import FloatArray
from .love_1d import complex_love_general, complex_love, effective_rigidity, effective_rigidity_general

ModeListType = Union[Tuple[np.ndarray], List[np.ndarray]]


@njit
def calc_tides(gravity: float, radius: float, density: float,
               shear_modulus: FloatArray, tidal_susceptibility: FloatArray,
               cmplx_compliance_bymode: ModeListType, heating_coeffs_bymode: ModeListType,
               torque_coeffs_bymode: ModeListType):
    """ Calculate Tidal Heating and Torques using the 2nd order Love number

    Parameters
    ----------
    gravity : float
        Surface gravity of layer [m s-2]
    radius : float
        Surface radius of layer [m]
    density : float
        Bulk density of layer [kg m-3]
    shear_modulus : FloatArray
        Rigidity of layer [Pa]
    tidal_susceptibility : FloatArray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]
    cmplx_compliance_bymode : ModeListType
        Complex Compliance calculated for each tidal mode (frequency), list of [Pa-1]
    heating_coeffs_bymode : ModeListType
        Tidal Heating coefficients calculated for each tidal mode
    torque_coeffs_bymode : ModeListType
        Tidal Torque coefficients calculated for each tidal mode

    Returns
    -------
    eff_rigid : np.ndarray
        Layer's effective rigidity
    tidal_heating : np.ndarray
        Layer's Tidal Heating [Watts]
    tidal_torque : np.ndarray
        Layer's Tidal zTorque [N m]
    tidal_heating_bymode : ModeListType
        List of layer's tidal heating for each tidal mode, list of [Watts]
    tidal_ztorque_bymode : ModeListType
        List of layer's tidal zTorque for each tidal mode, list of [N m]
    """

    # Calculate Effective Rigidity
    eff_rigid = effective_rigidity(shear_modulus, gravity, radius, density)

    # Calculate Complex Love Number
    cmplx_love_bymode = [complex_love(cmplx_comp, shear_modulus, eff_rigid)
                         for cmplx_comp in cmplx_compliance_bymode]

    # Calculate Tidal Heating
    tidal_heating_bymode = [-np.imag(love) * heating_coeff * tidal_susceptibility
                            for love, heating_coeff in zip(cmplx_love_bymode, heating_coeffs_bymode)]
    tidal_heating = tidal_heating_bymode[0].copy()
    for heating in tidal_heating_bymode[1:]:
        tidal_heating += heating

    # Calculate Tidal zTorque
    tidal_ztorque_bymode = [-np.imag(love) * torque_coeff * tidal_susceptibility
                            for love, torque_coeff in zip(cmplx_love_bymode, torque_coeffs_bymode)]
    tidal_torque = tidal_ztorque_bymode[0].copy()
    for torque in tidal_ztorque_bymode[1:]:
        tidal_torque += torque

    return eff_rigid, tidal_heating, tidal_torque, tidal_heating_bymode, tidal_ztorque_bymode


@njit
def calc_tides_general(gravity: float, radius: float, density: float,
               shear_modulus: np.ndarray, tidal_susceptibility: np.ndarray,
               cmplx_compliance_bymode: ModeListType, heating_coeffs_bymode: ModeListType,
               torque_coeffs_bymode: ModeListType, order_l: int = 2):
    """ Calculate Tidal Heating and Torques using the general Love number of order: order_l

    Parameters
    ----------
    gravity : float
        Surface gravity of layer [m s-2]
    radius : float
        Surface radius of layer [m]
    density : float
        Bulk density of layer [kg m-3]
    shear_modulus : np.ndarray
        Rigidity of layer [Pa]
    tidal_susceptibility : np.ndarray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]
    cmplx_compliance_bymode : ModeListType
        Complex Compliance calculated for each tidal mode (frequency), list of [Pa-1]
    heating_coeffs_bymode : ModeListType
        Tidal Heating coefficients calculated for each tidal mode
    torque_coeffs_bymode : ModeListType
        Tidal Torque coefficients calculated for each tidal mode
    order_l : int = 2
        First Fourier integer that sets the order for Love and Effective rigidity calculations

    Returns
    -------
    eff_rigid : np.ndarray
        Layer's effective rigidity
    tidal_heating : np.ndarray
        Layer's Tidal Heating [Watts]
    tidal_torque : np.ndarray
        Layer's Tidal zTorque [N m]
    tidal_heating_bymode : ModeListType
        List of layer's tidal heating for each tidal mode, list of [Watts]
    tidal_ztorque_bymode : ModeListType
        List of layer's tidal zTorque for each tidal mode, list of [N m]
    """

    # Calculate Effective Rigidity
    eff_rigid = effective_rigidity_general(shear_modulus, gravity, radius, density, order_l)

    # Calculate Complex Love Number
    cmplx_love_bymode = [complex_love_general(cmplx_comp, shear_modulus, eff_rigid, order_l)
                         for cmplx_comp in cmplx_compliance_bymode]

    # Calculate Tidal Heating
    tidal_heating_bymode = [-np.imag(love) * heating_coeff * tidal_susceptibility
                            for love, heating_coeff in zip(cmplx_love_bymode, heating_coeffs_bymode)]
    tidal_heating = tidal_heating_bymode[0].copy()
    for heating in tidal_heating_bymode[1:]:
        tidal_heating += heating

    # Calculate Tidal zTorque
    tidal_ztorque_bymode = [-np.imag(love) * torque_coeff * tidal_susceptibility
                            for love, torque_coeff in zip(cmplx_love_bymode, torque_coeffs_bymode)]
    tidal_torque = tidal_ztorque_bymode[0].copy()
    for torque in tidal_ztorque_bymode[1:]:
        tidal_torque += torque

    return eff_rigid, tidal_heating, tidal_torque, tidal_heating_bymode, tidal_ztorque_bymode