from typing import List, Tuple, Union

import numpy as np

from .love_1d import complex_love, complex_love_general, effective_rigidity, effective_rigidity_general
from ..performance import njit
from ..types import FloatArray
from ..dynamics import mode_types

from typing import Tuple

ModeListType = Union[Tuple[np.ndarray], List[np.ndarray]]

@njit
def collapse_terms(heating_terms: List[np.ndarray], torque_terms: List[np.ndarray], tidal_susceptibility: np.ndarray):
    """

    Parameters
    ----------
    heating_terms
    torque_terms
    tidal_susceptibility

    Returns
    -------

    """

    total_heating = heating_terms[0].copy()
    total_torque = torque_terms[0].copy()
    for heating_term, torque_term in zip(heating_terms[1:], torque_terms[1:]):
        total_heating += heating_term
        total_torque += torque_term

    total_heating *= tidal_susceptibility
    total_torque *= tidal_susceptibility

    return total_heating, total_torque


def calculate_tides(compliance: np.ndarray, viscosity: np.ndarray,
                    gravity: float, radius: float, density: float, tidal_susceptibility: FloatArray,
                    complex_compliance_func: object, complex_compliance_input: tuple,
                    semi_major_axis: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray,
                    orbital_frequency: np.ndarray, spin_frequency: np.ndarray,
                    use_nsr: bool = False, truncation: int = 2, order_l: int = 2):
    """ Calculate Tidal Heating and Torques using the 2nd order Love number

    Parameters
    ----------
    complex_compliance_func : object
        Complex compliance function to calculate love number at each mode.
    complex_compliance_input : tuple
        Other input used in the complex compliance calculation (other rheological parameters)
    compliance : FloatArray
        Compliance (inverse of rigidity) of the planet/layer's material [Pa-1]
    viscosity : FloatArray
        Viscosity of the planet/layer's material [Pa s]
    gravity : float
        Surface gravity of layer [m s-2]
    radius : float
        Surface radius of layer [m]
    density : float
        Bulk density of layer [kg m-3]
    tidal_susceptibility : FloatArray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]
    complex_compliance_func : object
        Complex compliance function to calculate love number at each mode.
    complex_compliance_input : tuple
        Other input used in the complex compliance calculation (other rheological parameters)
    eccentricity:
        Orbital Eccentricity
    inclination:
        Orbital Inclination (relative to the plane connecting satellite-to-host) [rads]
    orbital_frequency : FloatArray
        Orbital Mean Motion [rad s-1]
    spin_frequency : FloatArray
        Spin Frequency [rad s-1]
    use_nsr : bool = False
        Determines if non-synchronous rotation equations are used.
    truncation : int = 2
        Truncation to use on powers of e and I.
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

    if order_l < 2:
        raise ValueError
    elif order_l > 3:
        # TODO
        raise NotImplementedError

    # Calculate Effective Rigidity and Get Heating/Torque Equations
    shear_modulus = compliance**(-1)

    # Store heating and torque terms
    heating_terms = list()
    torque_terms = list()

    # cache complex compliances so repeat calculations are not made for the same mode.
    cached_comps = dict()

    for l in range(2, order_l):

        # Calculate the effective rigidity and multiplier at this l value
        eff_rigid = effective_rigidity_general(shear_modulus, gravity, radius, density, order_l=l)
        mode_calculator = mode_types[use_nsr][l][truncation]
        if l == 2:
            l_coeff = np.ones_like(semi_major_axis)
        else:
            # The -4 comes from that tidal_susceptibility already has a (R/a)^5 so l=3 should result in a coeff
            #  of (R/a)^2
            l_coeff = (radius / semi_major_axis)**(2*l - 4)

        mode_names, modes, freqs, heating_coeffs, ztorque_coeffs = \
            mode_calculator(orbital_frequency, spin_frequency, eccentricity, inclination)

        for mode_name, freq, heating_coeff, ztorque_coeff in zip(mode_names, freqs, heating_coeffs, ztorque_coeffs):

            # Calculate the complex compliance
            if mode_name not in cached_comps:
                cmplx_comp = complex_compliance_func(compliance, viscosity, freq, *complex_compliance_input)
                cached_comps[mode_name] = cmplx_comp
            else:
                cmplx_comp = cached_comps[mode_name]

            # Calculate the complex Love Number
            cmplx_love_number = complex_love_general(cmplx_comp, shear_modulus, eff_rigid, order_l=l)
            heating_terms.append(heating_coeff * -np.imag(cmplx_love_number) * l_coeff)
            torque_terms.append(ztorque_coeff * -np.imag(cmplx_love_number) * l_coeff)

    # Collapse all the tidal modes into a single number
    tidal_heating, tidal_torque = collapse_terms(heating_terms, torque_terms, tidal_susceptibility)

    return tidal_heating, tidal_torque
