from typing import Tuple, Dict, Union

import numpy as np

from .eccentricityFuncs import EccenOutput
from .inclinationFuncs import InclinOutput
from .love1d import effective_rigidity_general, complex_love_general
from .modeCalcHelper import eccentricity_functions_lookup, inclination_functions_lookup
from .universal_coeffs import get_universal_coeffs
from ..configurations import cache_numba
from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray, ComplexArray, NoneType

FreqSig = Tuple[int, int]
DissipTermsMix = Tuple[FloatArray, FloatArray, FloatArray, FloatArray]


@njit(cache=cache_numba)
def calculate_terms(spin_frequency: FloatArray, orbital_frequency: FloatArray,
                    semi_major_axis: FloatArray, radius: float,
                    eccentricity_results_byorderl: Dict[int, EccenOutput],
                    obliquity_results_byorderl: Dict[int, InclinOutput]):
    """ Calculate tidal dissipation terms and frequencies based on the current orbital and spin state.

    Requires the user to provide the correct eccentricity and inclination functions. This is generally done by the
    tides.calculate_terms wrapper.

    Parameters
    ----------
    orbital_frequency : FloatArray
        Orbit's orbital frequency around tidal host [rad s-1]
    spin_frequency : FloatArray
        World's rotation frequency [rad s-1]
    semi_major_axis : FloatArray
        Orbital separation [m]
    radius : float
        Planet's radius (used to calculate R/a) [m]
    eccentricity_results_byorderl : Dict[int, EccenOutput]
        Pre-calculated eccentricity results using the .modeCalcHelper functions
        Stored as a numba-safe dict of order-l: eccen_result.
    obliquity_results_byorderl : Dict[int, EccenOutput]
        Pre-calculated obliquity results using the .modeCalcHelper functions
        Stored as a numba-safe dict of order-l: obliquity_result.

    Returns
    -------
    unique_frequencies : Dict[FreqSig, FloatArray]
        Unique frequency signature and calculate frequency [rad s-1]
        Frequencies are integer combinations of orbital motion and spin-rate.
    results_by_frequency : Dict[FreqSig, Dict[int, DissipTermsMix]]
        Tidal dissipation terms stored for each unique frequency and tidal harmonic order-l.
        At each frequency dissipation terms are stored as:
            heating_term, dUdM_term, dUdw_term, dUdO_term
            These still need to be multiplied by the -Im[k_l] and tidal susceptibility.

    See Also
    --------
    TidalPy.tides.dissipation.calculate_terms
    """

    # Determine the maximum order l from the eccentricity and obliquity results
    max_order_l_eccen = 1 + len(eccentricity_results_byorderl)
    max_order_l_obliquity = 1 + len(obliquity_results_byorderl)
    if max_order_l_eccen != max_order_l_obliquity:
        # The maximum order l should be the same for obliquity and eccentricity results.
        #    Mismatch in eccentricity obliquity function?
        raise Exception
    max_order_l = max_order_l_eccen

    # Storage for results by unique frequency signature. Must provide fake data structures and values so that numba
    #    has an idea on how to compile.
    eccen_res_for_fake = eccentricity_results_byorderl[2][0][0]
    obliq_res_for_fake = obliquity_results_byorderl[2][(0, 1)]
    fake_result = orbital_frequency * spin_frequency * eccen_res_for_fake * obliq_res_for_fake
    results_by_frequency = {
        (-100, -100): {
            0: (fake_result, fake_result, fake_result, fake_result)
        }
    }
    unique_frequencies = {
        (-100, -100): orbital_frequency - spin_frequency
    }

    # Distance Scale used in the distance multiplier
    distance_scale = radius / semi_major_axis

    # Main calculate loop
    for order_l in range(2, max_order_l + 1):

        # Get universal coefficients for this order l
        universal_coeff_by_m = get_universal_coeffs(order_l)

        # Eccentricity and obliquity results are pre-calculated for each tidal order l.
        #    Get result for this order l
        obliquity_results = obliquity_results_byorderl[order_l]
        eccentricity_results = eccentricity_results_byorderl[order_l]

        # Order l controls the distance multiplier. The minus 4 comes from the tidal susceptibility already carrying
        #    (R / a)^5 so (R / a)^(2l + 1) --> (R / a)^(2l - 4)
        distance_multiplier = distance_scale**(2 * order_l - 4)

        for (m, p), obliquity_terms in obliquity_results.items():
            # Pull out the m and p integers from the non-zero obliquity terms

            # Pull out universal coefficient
            #    The Tidal Susceptibility carries a factor of (3 / 2) already. So we need to divide the uni_coeff
            #    by that much to ensure it is not double counted.
            uni_coeff = distance_multiplier * universal_coeff_by_m[m] / 1.5

            for q, eccentricity_terms in eccentricity_results[p].items():
                # Pull out q integer from the non-zero eccentricity terms

                # Multiply eccentricity terms by the obliquity terms
                eccen_obliq = obliquity_terms * eccentricity_terms
                uni_multiplier = uni_coeff * eccen_obliq

                # Calculate tidal mode, frequency, and sign
                n_coeff = (order_l - 2 * p + q)

                # Skip terms that we know are zero (0 * n - 0 * m will always be 0 freq and -imk(0) == 0)
                if n_coeff == 0 and m == 0:
                    continue

                mode = n_coeff * orbital_frequency - m * spin_frequency
                mode_sign = np.sign(mode)
                mode_frequency = np.abs(mode)

                # Determine the frequency signature used to store unique frequencies
                freq_sig = (n_coeff, -m)

                # Calculate coefficients for heating and potential derivatives
                heating_term = uni_multiplier * mode_frequency
                dUdM_term = uni_multiplier * n_coeff * mode_sign
                dUdw_term = uni_multiplier * (order_l - 2 * p) * mode_sign
                dUdO_term = uni_multiplier * m * mode_sign

                # The tidal heating and potential derivatives should also be multiplied by the love number calculated
                #    for that mode. But, the Love number only cares about the complex compliance which in turn only requires
                #    frequencies (does not care about modes or order-l). Therefore, we can store everything in terms of
                #    *unique* frequencies.
                if freq_sig not in unique_frequencies:
                    # New unique frequency
                    unique_frequencies[freq_sig] = mode_frequency
                    results_by_frequency[freq_sig] = {order_l: (heating_term, dUdM_term, dUdw_term, dUdO_term)}
                else:
                    if order_l in results_by_frequency[freq_sig]:
                        # Previous results found at this frequency and this order_l.
                        #     Add new results to previous one.
                        heating_term_old, dUdM_term_old, dUdw_term_old, dUdO_term_old = \
                            results_by_frequency[freq_sig][order_l]
                        heating_term_new = heating_term_old + heating_term
                        dUdM_term_new = dUdM_term_old + dUdM_term
                        dUdw_term_new = dUdw_term_old + dUdw_term
                        dUdO_term_new = dUdO_term_old + dUdO_term

                        # Replace old terms with the new ones
                        results_by_frequency[freq_sig][order_l] = \
                            (heating_term_new, dUdM_term_new, dUdw_term_new, dUdO_term_new)
                    else:
                        # No results at this frequency for this order-l
                        results_by_frequency[freq_sig][order_l] = (heating_term, dUdM_term, dUdw_term, dUdO_term)

    # Delete those fake results
    del results_by_frequency[(-100, -100)]
    del unique_frequencies[(-100, -100)]

    return unique_frequencies, results_by_frequency

@njit(cache=cache_numba)
def collapse_modes(gravity: float, radius: float, density: float,
                   shear_modulus: Union[NoneType, FloatArray], tidal_scale: float,
                   tidal_host_mass: float,
                   tidal_susceptibility: FloatArray,
                   complex_compliance_by_frequency: Tuple[ComplexArray, ...],
                   tidal_terms_by_frequency: Dict[FreqSig, Dict[int, DissipTermsMix]],
                   max_order_l: int, cpl_ctl_method: bool = False) -> \
        Tuple[FloatArray, FloatArray, FloatArray, FloatArray, ComplexArray, FloatArray]:
    """ Collapses the tidal terms calculated by calculate_modes() combined with rheological information from the layer.

    Parameters
    ----------
    gravity : float
        Acceleration due to gravity at the surface of the layer or planet [m s-2]
    radius : float
        Surface radius of the layer or planet [m]
    density : float
        Bulk density of the layer or planet [kg m-3]
    shear_modulus : FloatArray
        Effective shear modulus of the layer or planet [Pa]
    tidal_scale : float
        Multiplier that ranges from 0 to 1 that scales the imaginary portion of the complex love number to simulate
            only a portion of a planet's volume is contributing to the tidal dissipation.
    tidal_host_mass : float
        Mass of the tidal host [kg]
        Needed to offset the tidal susceptibility for the potential derivatives.
    tidal_susceptibility : FloatArray
        Tidal susceptibility, defined as (3/2) G M_host^2 R_target^5 / a^6 [N m]
    complex_compliance_by_frequency : Tuple[ComplexArray]
        The complex compliance for the layer or planet calculated at each unique tidal frequency [Pa-1]
    tidal_terms_by_frequency : Dict[FreqSig, Dict[int, DissipTerms]]
        Each dissipation term: E^dot, dUdM, dUdw, dUdO; calculated for each unique tidal frequency and order-l
    max_order_l : int
        Max tidal order level.
    cpl_ctl_method : bool = False
        Changes functionality based on if the method is ctl cpl or neither
        See tides.SimpleTides.mode_collapse

    Returns
    -------
    tidal_heating : FloatArray
        Tidal heating [W]
        This could potentially restricted to a layer or for an entire planet.
    dUdM : FloatArray
        Tidal potential derivative with respect to the mean anomaly [J kg-1 radians-1]
        This could potentially restricted to a layer or for an entire planet.
    dUdw : FloatArray
        Tidal potential derivative with respect to the pericentre [J kg-1 radians-1]
        This could potentially restricted to a layer or for an entire planet.
    dUdO : FloatArray
        Tidal potential derivative with respect to the planet's node [J kg-1 radians-1]
        This could potentially restricted to a layer or for an entire planet.
    love_number : ComplexArray
        The complex Love number [unitless]
        This could potentially restricted to a layer or for an entire planet.
    negative_imk : FloatArray
        The negative of the imaginary part of the complex Love number, scaled by the planet's tidal scale
        This could potentially restricted to a layer or for an entire planet.

    See Also
    --------
    TidalPy.tides.dissipation.calculate_modes
    TidalPy.tides.tides.Tides

    """

    # Check that the number of tidal terms (by freq) should be equal to number of complex compliances (by freq)
    if len(tidal_terms_by_frequency) != len(complex_compliance_by_frequency):
        # Mismatch between unique tidal modes and tidal terms.
        #    Ensure you are using the complex compliance from the correct object (different planets may have different
        #    number of modes)
        raise Exception

    tidal_heating_terms = list()
    dUdM_terms = list()
    dUdw_terms = list()
    dUdO_terms = list()
    love_number_terms = list()
    negative_imk_terms = list()
    signatures = list()

    for tidal_order_l in range(2, max_order_l + 1):

        # Calculate the effective rigidity which does not use the complex compliance.
        # TODO: Should this use the static shear modulus or the effective (post-melting/heating) shear?
        #    Right now it is using the effective rigidity. If it only needs the static then this could be done outside
        #    this function.
        if cpl_ctl_method:
            effective_rigidity = effective_rigidity_general(shear_modulus, gravity, radius, density,
                                                            order_l=2)
        else:
            if shear_modulus is None:
                # Numba does not support specific exception classes.
                raise Exception('Shear modulus is required for non-CTL/CPL models. Set it to a fake float e.g., 1.')

            effective_rigidity = effective_rigidity_general(shear_modulus, gravity, radius, density,
                                                            order_l=tidal_order_l)

        # Pull out the already computed tidal heating and potential terms each frequency
        freq_i = 0
        for unique_freq_signature, tidal_terms in tidal_terms_by_frequency.items():

            # Many higher-order frequencies do not have lower-order l results, so skip if this order-l is not in
            #    this tidal terms
            if tidal_order_l in tidal_terms:
                # Otherwise: Pull out the complex compliance at this frequency
                complex_compliance = complex_compliance_by_frequency[freq_i]

                if cpl_ctl_method:
                    # In the CTL/CPL method, the complex love number is passed in as the complex compliance
                    complex_love = complex_compliance
                else:
                    # Calculate the complex Love number, which also changes functional form for each order_l
                    complex_love = complex_love_general(complex_compliance, shear_modulus, effective_rigidity,
                                                        order_l=tidal_order_l)

                # Scale the Love number by the layer's contribution
                # TODO: Should the tidal scale affect the Re(k) as well as the Im(k)?
                complex_love = np.real(complex_love) + (1.0j * np.imag(complex_love) * tidal_scale)
                neg_imk = -np.imag(complex_love)

                # The tidal potential carries one fewer dependence on the tidal host mass that is built into the
                #    tidal susceptibility. Divide that out now.
                neg_imk_potential = neg_imk / tidal_host_mass

                # Pull out the tidal terms pre-calculated for this unique frequency. See Tides.update_orbit_spin
                heating_term, dUdM_term, dUdw_term, dUdO_term = tidal_terms[tidal_order_l]

                # Store results
                tidal_heating_terms.append(heating_term * neg_imk)
                dUdM_terms.append(dUdM_term * neg_imk_potential)
                dUdw_terms.append(dUdw_term * neg_imk_potential)
                dUdO_terms.append(dUdO_term * neg_imk_potential)
                love_number_terms.append(complex_love)
                negative_imk_terms.append(neg_imk)
                signatures.append(unique_freq_signature)

            # Increment the frequency index rather or not the tidal order l was found for this tidal terms.
            freq_i += 1

    # Collapse Modes
    # FIXME: Njit did not like sum( ), so doing separate loop for these for now...
    tidal_heating = tidal_heating_terms[0]
    dUdM = dUdM_terms[0]
    dUdw = dUdw_terms[0]
    dUdO = dUdO_terms[0]
    love_number = love_number_terms[0]
    negative_imk = negative_imk_terms[0]

    for term_i in range(1, len(tidal_heating_terms)):

        tidal_heating += tidal_heating_terms[term_i]
        dUdM += dUdM_terms[term_i]
        dUdw += dUdw_terms[term_i]
        dUdO += dUdO_terms[term_i]
        love_number += love_number_terms[term_i]
        negative_imk += negative_imk_terms[term_i]

    tidal_heating *= tidal_susceptibility
    dUdM *= tidal_susceptibility
    dUdw *= tidal_susceptibility
    dUdO *= tidal_susceptibility

    return tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk


def find_mode_manipulators(max_order_l: int = 2, eccentricity_truncation_lvl: int = 8, use_obliquity: bool = True):

    # Find Eccentricity and Obliquity Functions
    eccentricity_func = eccentricity_functions_lookup[eccentricity_truncation_lvl][max_order_l]
    inclination_func = inclination_functions_lookup[use_obliquity][max_order_l]

    return calculate_terms, collapse_modes, eccentricity_func, inclination_func