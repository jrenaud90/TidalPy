from typing import List, Tuple, Dict

import numpy as np

from .love1d import complex_love_general, effective_rigidity_general
from ..performance import njit
from ..types import FloatArray
from ..constants import G


@njit
def calc_tidal_susceptibility_reduced(host_mass: float, target_radius: float) -> np.ndarray:
    """ Calculate the reduced tidal susceptibility for a given target radius, host mass, and their separation.

    Parameters
    ----------
    host_mass : float
        Mass of central host [kg]
    target_radius : float
        Radius of target body [m]

    Returns
    -------
    tidal_susceptibility_reduced : np.ndarray
        Reduced Tidal Susceptibility [kg m8 s-2]
    """

    tidal_susceptibility_reduced = (3. / 2.) * G * host_mass**2 * target_radius**5

    return tidal_susceptibility_reduced





@njit
def _order_l_coeffs(max_order_l: int, semi_major_axis: np.ndarray, world_radius: float):

    coeffs = list()
    for l in range(2, max_order_l+1):
        if l == 2:
            coeff = np.ones_like(semi_major_axis)
        else:
            # The -4 comes from that tidal_susceptibility already has a (R/a)^5 so l=3 should result in a coeff
            #  of (R/a)^2
            coeff = (world_radius / semi_major_axis)**(2 * l - 4)

        coeffs.append(coeff)

    return coeffs


@njit
def _collapse_terms(heating_terms: List[np.ndarray], torque_terms: List[np.ndarray],
                    tidal_susceptibility: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """ Collapses multiple tidal mode results down into a signle value.

    This is separated out from calculate_tides so it can be wrapped in njit

    Parameters
    ----------
    heating_terms : List[np.ndarray]
        List of tidal heating terms (each should be unitless as tidal_susceptibility carries all units)
    torque_terms : List[np.ndarray]
        List of tidal torque terms (each should be unitless as tidal_susceptibility carries all units)
    tidal_susceptibility : np.ndarray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]

    Returns
    -------
    total_heating : np.ndarray
        Total Tidal Heating [W]
    total_torque : np.ndarray
        Total Tidal Torque [N m]
    """

    total_heating = np.copy(heating_terms[0])
    total_ztorque = np.copy(torque_terms[0])
    for heating_term, ztorque_term in zip(heating_terms[1:], torque_terms[1:]):
        total_heating += heating_term
        total_ztorque += ztorque_term

    total_heating *= tidal_susceptibility
    total_ztorque *= tidal_susceptibility

    return total_heating, total_ztorque


def calculate_tides(shear_modulus: np.ndarray, viscosity: np.ndarray,
                    gravity: float, radius: float, density: float, tidal_susceptibility: FloatArray,
                    complex_compliance_func: object, complex_compliance_input: tuple,
                    semi_major_axis: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray,
                    orbital_frequency: np.ndarray, spin_frequency: np.ndarray,
                    tidal_volume_fraction: float = 1.,
                    use_nsr: bool = False, truncation: int = 2, order_l: int = 2) -> \
        Tuple[np.ndarray, np.ndarray, List[np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """ Calculate Tidal Heating and Torques using a generalized complex Love number.

    Parameters
    ----------
    shear_modulus : np.ndarray
        Shear modulus of the planet/layer's material [Pa]
    viscosity : np.ndarray
        Viscosity of the planet/layer's material [Pa s]
    gravity : float
        Surface gravity of layer [m s-2]
    radius : float
        Surface radius of layer [m]
    density : float
        Bulk density of layer [kg m-3]
    tidal_susceptibility : np.ndarray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]
    complex_compliance_func : object
        Complex compliance function to calculate love number at each mode.
    complex_compliance_input : tuple
        Other input used in the complex compliance calculation (other rheological parameters)
    semi_major_axis : np.ndarray
        Semi-major axis [m]
    eccentricity : np.ndarray
        Orbital Eccentricity
    inclination : np.ndarray
        Orbital Inclination (relative to the plane connecting satellite-to-host) [rads]
    orbital_frequency : FloatArray
        Orbital Mean Motion [rad s-1]
    spin_frequency : FloatArray
        Spin Frequency [rad s-1]
    tidal_volume_fraction : float
        Tidal volume fraction = tidal volume / total planet volume
        Way to estimate multi-layer dissipation
    use_nsr : bool = False
        Determines if non-synchronous rotation equations are used.
    truncation : int = 2
        Truncation to use on powers of e and I.
    order_l : int = 2
        First Fourier integer that sets the order for Love and Effective rigidity calculations

    Returns
    -------
    tidal_heating : np.ndarray
        Layer's Tidal Heating [W]
    tidal_torque : np.ndarray
        Layer's Tidal zTorque [N m]
    love_numbers : List[np.ndarray, ...]
        The complex Love numbers for each tidal mode
    tidal_frequencies : List[np.ndarray, ...]
        Tidal frequencies [rads s-1]
    cached_complex_comps : List[np.ndarray, ...]
        Complex compliances for each frequency [Pa-1]
    """

    if order_l < 2:
        raise ValueError
    elif order_l > 3:
        # TODO
        raise NotImplementedError

    # Calculate Effective Rigidity and Get Heating/Torque Equations.
    compliance = 1. / shear_modulus

    # Store heating and torque terms as well as the love numbers.
    heating_terms = list()
    ztorque_terms = list()
    love_numbers = list()

    # Cache complex compliances so repeat calculations are not made for the same mode.
    cached_complex_comps = dict()
    tidal_frequencies = dict()

    for l in range(2, order_l+1):

        # Calculate the effective rigidity and multiplier at this l value
        eff_rigid = effective_rigidity_general(shear_modulus, gravity, radius, density, order_l=l)
        mode_calculator = mode_types[use_nsr][l][truncation]
        if l == 2:
            l_coeff = np.ones_like(semi_major_axis)
        else:
            # The -4 comes from that tidal_susceptibility already has a (R/a)^5 so l=3 should result in a coeff
            #  of (R/a)^2
            l_coeff = (radius / semi_major_axis)**(2*l - 4)

        # Scale by tvf
        l_coeff *= tidal_volume_fraction

        mode_names, modes, freqs, heating_subterms, ztorque_subterms = \
            mode_calculator(orbital_frequency, spin_frequency, eccentricity, inclination)

        for mode_name, freq, heating_coeff, ztorque_coeff in zip(mode_names, freqs, heating_subterms, ztorque_subterms):

            # Calculate the complex compliance
            if mode_name not in cached_complex_comps:
                # This particular frequency has not been calculated before.
                cmplx_comp = complex_compliance_func(compliance, viscosity, freq, *complex_compliance_input)
                cached_complex_comps[mode_name] = cmplx_comp
                tidal_frequencies[mode_name] = freq
            else:
                cmplx_comp = cached_complex_comps[mode_name]

            # Calculate the complex Love Number
            cmplx_love_number = complex_love_general(cmplx_comp, shear_modulus, eff_rigid, order_l=l)
            neg_imk = -np.imag(cmplx_love_number)

            love_numbers.append(cmplx_love_number)
            heating_terms.append(heating_coeff * neg_imk * l_coeff)
            ztorque_terms.append(ztorque_coeff * neg_imk * l_coeff)

    # Collapse all the tidal modes into a single number
    tidal_heating, tidal_ztorque = _collapse_terms(heating_terms, ztorque_terms, tidal_susceptibility)

    return tidal_heating, tidal_ztorque, love_numbers, tidal_frequencies, cached_complex_comps


def calculate_tides_withmodes(shear_modulus: np.ndarray, viscosity: np.ndarray,
                              gravity: float, radius: float, density: float, tidal_susceptibility: FloatArray,
                              complex_compliance_func: object, complex_compliance_input: tuple,
                              semi_major_axis: np.ndarray,
                              mode_names_byorderl: List[List[str]], modes_byorderl: List[List[np.ndarray]],
                              freqs_byorderl: List[List[np.ndarray]], heating_subterms_byorderl: List[List[np.ndarray]],
                              ztorque_subterms_byorderl: List[List[np.ndarray]],
                              tidal_volume_fraction: float = 1., order_l: int = 2) -> \
        Tuple[np.ndarray, np.ndarray, List[np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """ Calculate Tidal Heating and Torques using a generalized complex Love number and user provided mode information.

    Parameters
    ----------
    shear_modulus : np.ndarray
        Shear modulus of the planet/layer's material [Pa]
    viscosity : np.ndarray
        Viscosity of the planet/layer's material [Pa s]
    gravity : float
        Surface gravity of layer [m s-2]
    radius : float
        Surface radius of layer [m]
    density : float
        Bulk density of layer [kg m-3]
    tidal_susceptibility : np.ndarray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]
    complex_compliance_func : object
        Complex compliance function to calculate love number at each mode.
    complex_compliance_input : tuple
        Other input used in the complex compliance calculation (other rheological parameters)
    semi_major_axis : np.ndarray
        Semi-major axis [m]
    mode_names_byorderl : List[List[str]]
        List of easy to read names  for each mode, and over each order-l
            "n" == orbital motion frequency
            "o" == spin frequency
    modes_byorderl : List[List[np.ndarray]]
        List of calculated modes for each mode, and over each order-l [rads s-1]
    freqs_byorderl : List[List[np.ndarray]]
        List of calculated frequencies (absolute value of modes) for each mode, and over each order-l [rads s-1]
    heating_subterms_byorderl: List[List[np.ndarray]]
        List of tidal heating coefficients for each mode, and over each order-l [N m]
    ztorque_subterms_byorderl: List[List[np.ndarray]]
        List of tidal torque coefficients for each mode, and over each order-l [N m]
    tidal_volume_fraction : float
        Tidal volume fraction = tidal volume / total planet volume
        Way to estimate multi-layer dissipation
    order_l : int = 2
        First Fourier integer that sets the order for Love and Effective rigidity calculations

    Returns
    -------
    tidal_heating : np.ndarray
        Layer's Tidal Heating [W]
    tidal_torque : np.ndarray
        Layer's Tidal zTorque [N m]
    love_numbers : List[np.ndarray, ...]
        The complex Love numbers for each tidal mode
    tidal_frequencies : List[np.ndarray, ...]
        Tidal frequencies [rads s-1]
    cached_complex_comps : List[np.ndarray, ...]
        Complex compliances for each frequency [Pa-1]
    """

    if order_l < 2:
        raise ValueError
    elif order_l > 3:
        # TODO
        raise NotImplementedError

    # Calculate Effective Rigidity and Get Heating/Torque Equations.
    compliance = 1. / shear_modulus

    # Store heating and torque terms as well as the love numbers.
    heating_terms = list()
    ztorque_terms = list()
    love_numbers = list()

    # Cache complex compliances so repeat calculations are not made for the same mode.
    cached_complex_comps = dict()
    tidal_frequencies = dict()

    for l in range(2, order_l+1):

        # Calculate the effective rigidity and multiplier at this l value
        eff_rigid = effective_rigidity_general(shear_modulus, gravity, radius, density, order_l=l)
        if l == 2:
            l_coeff = np.ones_like(semi_major_axis)
        else:
            # The -4 comes from that tidal_susceptibility already has a (R/a)^5 so l=3 should result in a coeff
            #  of (R/a)^2
            l_coeff = (radius / semi_major_axis)**(2*l - 4)

        # Scale by tvf
        l_coeff *= tidal_volume_fraction

        mode_names, modes, freqs, heating_subterms, ztorque_subterms = \
            mode_names_byorderl[order_l], modes_byorderl[order_l], freqs_byorderl[order_l], \
            heating_subterms_byorderl[order_l], ztorque_subterms_byorderl[order_l]

        for mode_name, freq, heating_coeff, ztorque_coeff in zip(mode_names, freqs, heating_subterms, ztorque_subterms):

            # Calculate the complex compliance
            if mode_name not in cached_complex_comps:
                # This particular frequency has not been calculated before.
                cmplx_comp = complex_compliance_func(compliance, viscosity, freq, *complex_compliance_input)
                cached_complex_comps[mode_name] = cmplx_comp
                tidal_frequencies[mode_name] = freq
            else:
                cmplx_comp = cached_complex_comps[mode_name]

            # Calculate the complex Love Number
            cmplx_love_number = complex_love_general(cmplx_comp, shear_modulus, eff_rigid, order_l=l)
            neg_imk = -np.imag(cmplx_love_number)

            love_numbers.append(cmplx_love_number)
            heating_terms.append(heating_coeff * neg_imk * l_coeff)
            ztorque_terms.append(ztorque_coeff * neg_imk * l_coeff)

    # Collapse all the tidal modes into a single number
    tidal_heating, tidal_ztorque = _collapse_terms(heating_terms, ztorque_terms, tidal_susceptibility)

    return tidal_heating, tidal_ztorque, love_numbers, tidal_frequencies, cached_complex_comps


def calculate_tides_withlove(shear_modulus: np.ndarray, viscosity: np.ndarray,
                             gravity: float, world_radius: float, density: float, tidal_susceptibility: FloatArray,
                             complex_love_numbers_byorderl: List[Dict[str, np.ndarray]],
                             semi_major_axis: np.ndarray,
                             mode_names_byorderl: List[List[str]],
                             heating_subterms_byorderl: List[List[np.ndarray]],
                             ztorque_subterms_byorderl: List[List[np.ndarray]],
                             tidal_volume_fraction: float = 1., order_l: int = 2) -> \
        Tuple[np.ndarray, np.ndarray, List[np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """ Calculate Tidal Heating and Torques using a generalized complex Love number and user provided mode information.

    Parameters
    ----------
    shear_modulus : np.ndarray
        Shear modulus of the planet/layer's material [Pa]
    viscosity : np.ndarray
        Viscosity of the planet/layer's material [Pa s]
    gravity : float
        Surface gravity of layer [m s-2]
    radius : float
        Surface radius of layer [m]
    density : float
        Bulk density of layer [kg m-3]
    tidal_susceptibility : np.ndarray
        Planet's tidal susceptibility (3/2 * G * M_host^2 * R_planet^5 / a^6) [N m]
    complex_compliance_func : object
        Complex compliance function to calculate love number at each mode.
    complex_compliance_input : tuple
        Other input used in the complex compliance calculation (other rheological parameters)
    semi_major_axis : np.ndarray
        Semi-major axis [m]
    eccentricity : np.ndarray
        Orbital Eccentricity
    inclination : np.ndarray
        Orbital Inclination (relative to the plane connecting satellite-to-host) [rads]
    mode_names_byorderl : List[List[str]]
        List of easy to read names  for each mode, and over each order-l
            "n" == orbital motion frequency
            "o" == spin frequency
    modes_byorderl : List[List[np.ndarray]]
        List of calculated modes for each mode, and over each order-l [rads s-1]
    freqs_byorderl : List[List[np.ndarray]]
        List of calculated frequencies (absolute value of modes) for each mode, and over each order-l [rads s-1]
    heating_subterms_byorderl: List[List[np.ndarray]]
        List of tidal heating coefficients for each mode, and over each order-l [N m]
    ztorque_subterms_byorderl: List[List[np.ndarray]]
        List of tidal torque coefficients for each mode, and over each order-l [N m]
    tidal_volume_fraction : float
        Tidal volume fraction = tidal volume / total planet volume
        Way to estimate multi-layer dissipation
    use_nsr : bool = False
        Determines if non-synchronous rotation equations are used.
    truncation : int = 2
        Truncation to use on powers of e and I.
    order_l : int = 2
        First Fourier integer that sets the order for Love and Effective rigidity calculations

    Returns
    -------
    tidal_heating : np.ndarray
        Layer's Tidal Heating [W]
    tidal_torque : np.ndarray
        Layer's Tidal zTorque [N m]
    love_numbers : List[np.ndarray, ...]
        The complex Love numbers for each tidal mode
    tidal_frequencies : List[np.ndarray, ...]
        Tidal frequencies [rads s-1]
    cached_complex_comps : List[np.ndarray, ...]
        Complex compliances for each frequency [Pa-1]
    """

    if order_l < 2:
        raise ValueError
    elif order_l > 3:
        # TODO
        raise NotImplementedError

    # Store heating and torque terms as well as the love numbers.
    heating_terms = list()
    ztorque_terms = list()

    for l in range(2, order_l+1):

        heating_subterms_bymode = heating_subterms_byorderl[l]
        ztorque_subterms_bymode = ztorque_subterms_byorderl[l]
        love_numbers_bymode = complex_love_numbers_byorderl[l]
        mode_names = mode_names_byorderl[l]

        # Calculate the multiplier at this l value
        if l == 2:
            l_coeff = np.ones_like(semi_major_axis)
        else:
            # The -4 comes from that tidal_susceptibility already has a (R/a)^5 so l=3 should result in a coeff
            #  of (R/a)^2
            l_coeff = (world_radius / semi_major_axis)**(2*l - 4)

        # Scale by tvf
        l_coeff *= tidal_volume_fraction

        for mode_name, heating_subterm, ztorque_subterm in \
                zip(mode_names, heating_subterms_bymode, ztorque_subterms_bymode):

            # Calculate the complex Love Number
            cmplx_love_number = love_numbers_bymode[mode_name]
            neg_imk = -np.imag(cmplx_love_number)

            heating_terms.append(heating_subterm * neg_imk * l_coeff)
            ztorque_terms.append(ztorque_subterm * neg_imk * l_coeff)

    # Collapse all the tidal modes into a single number
    tidal_heating, tidal_ztorque = _collapse_terms(heating_terms, ztorque_terms, tidal_susceptibility)

    return tidal_heating, tidal_ztorque, love_numbers, tidal_frequencies, cached_complex_comps