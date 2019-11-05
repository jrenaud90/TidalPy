from typing import List, Tuple, Union

import numpy as np

from .love_1d import complex_love, complex_love_general, effective_rigidity, effective_rigidity_general
from ..performance import njit
from ..types import FloatArray
from ..dynamics import mode_types


ModeListType = Union[Tuple[np.ndarray], List[np.ndarray]]


@njit
def calculate_tides(complex_compliance_func: object, complex_compliance_input: tuple,
                    compliance: FloatArray, viscosity: FloatArray,
                    gravity: float, radius: float, density: float, tidal_susceptibility: FloatArray,
                    eccentricity: FloatArray, inclination: FloatArray,
                    orbital_frequency: FloatArray, spin_frequency: FloatArray = None,
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
        raise NotImplementedError

    # Calculate Effective Rigidity and Get Heating/Torque Equations
    shear_modulus = compliance**(-1)
    eff_rigid_2 = effective_rigidity(shear_modulus, gravity, radius, density)
    calc_mode_2 = mode_types[use_nsr][2][truncation]

    eff_rigid_3 = None
    calc_mode_3 = None
    if order_l == 3:
        raise NotImplementedError # FIXME
        eff_rigid_3 = effective_rigidity_general(shear_modulus, gravity, radius, density, order_l=3)
        calc_mode_3 = mode_types[use_nsr][3][truncation]

    # Calculate Tidal Modes for l=2
    modes, freqs, heating_coeffs, ztorque_coeffs = \
        calc_mode_2(orbital_frequency, spin_frequency, eccentricity, inclination)

    calculated_freqs = list()
    calculated_comps = list()
    heating_terms = list()
    torque_terms = list()
    for freq, heating_coeff, torque_coeff in zip(freqs, heating_coeffs, ztorque_coeffs):

        # Store calculated frequencies so they can be checked for other l != 2 (the compliances will be the same)
        calculated_freqs.append(freq)
        cmplx_comp = complex_compliance_func(compliance, viscosity, freq, *complex_compliance_input)
        calculated_comps.append(cmplx_comp)

        # Calculate Love Number, then heating and torque
        cmplx_love_number = complex_love(cmplx_comp, shear_modulus, eff_rigid_2)
        heating_terms.append(heating_coeff * -np.imag(cmplx_love_number))
        torque_terms.append(torque_coeff * -np.imag(cmplx_love_number))

    if calc_mode_3 is not None:
        modes_3, freqs_3, heating_coeffs_3, ztorque_coeffs_3 = \
            calc_mode_3(orbital_frequency, spin_frequency, eccentricity, inclination)
        for freq_3, heating_coeff_3, torque_coeff_3 in zip(freqs_3, heating_coeffs_3, ztorque_coeffs_3):

            if freq_3 in calculated_freqs:
                cmplx_comp = calculated_comps[calculated_freqs.index(freq_3)]
            else:
                calculated_freqs.append(freq_3)
                cmplx_comp = complex_compliance_func(compliance, viscosity, freq_3, *complex_compliance_input)
                calculated_comps.append(cmplx_comp)

            # Calculate Love Number, then heating and torque
            order_3_scale = (radius / semi_major_axis)**2
            cmplx_love_number = complex_love_general(cmplx_comp, shear_modulus, eff_rigid_3, order_l=3)
            heating_terms.append(heating_coeff * -np.imag(cmplx_love_number) * order_3_scale)
            torque_terms.append(torque_coeff * -np.imag(cmplx_love_number) * order_3_scale)


    total_heating = heating_terms[0].copy()
    total_torque = torque_terms[0].copy()
    for heating_term, torque_term in zip(heating_terms[1:], torque_terms[1:]):
        total_heating += heating_term
        total_torque += torque_term

    total_heating *= tidal_susceptibility
    total_torque *= tidal_susceptibility

    return tidal_heating, tidal_torque, heating_terms, torque_terms
