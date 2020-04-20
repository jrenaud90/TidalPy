from typing import Tuple, Dict, List

import numpy as np
from numba import njit
from scipy.constants import G

from .eccentricityFuncs import eccentricity_truncations
from .inclinationFuncs import inclination_functions
from .love1d import effective_rigidity_general, complex_love_general
from .universal_coeffs import get_universal_coeffs
from ..types import FloatArray, ComplexArray

FreqSig = Tuple[int, int]
DissipTermsFloat = Tuple[float, float, float, float]
DissipTermsArray = Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
DissipTermsMix = Tuple[FloatArray, FloatArray, FloatArray, FloatArray]

@njit
def calc_tidal_susceptibility(host_mass: float, target_radius: float, semi_major_axis: FloatArray) -> FloatArray:
    """ Calculate the tidal susceptibility for a given target radius, host mass, and their separation.

    Parameters
    ----------
    host_mass : float
        Mass of central host [kg]
    target_radius : float
        Radius of target body [m]
    semi_major_axis : FloatArray
        Semi-major axis [m]

    Returns
    -------
    tidal_susceptibility : FloatArray
        Tidal Susceptibility [N m]
    """

    tidal_susceptibility = (3. / 2.) * G * host_mass**2 * target_radius**5 / semi_major_axis**6

    return tidal_susceptibility

@njit
def calc_tidal_susceptibility_reduced(host_mass: float, target_radius: float) -> float:
    """ Calculate the tidal susceptibility (reduced) for a given target radius and host mass.

    The reduced tidal susceptibility excludes the semi-major axis.

    Parameters
    ----------
    host_mass : float
        Mass of central host [kg]
    target_radius : float
        Radius of target body [m]

    Returns
    -------
    tidal_susceptibility_reduced : np.ndarray
        Tidal Susceptibility [N m]
    """

    tidal_susceptibility_reduced = (3. / 2.) * G * host_mass**2 * target_radius**5

    return tidal_susceptibility_reduced

def calculate_terms(orbital_frequency: FloatArray, spin_frequency: FloatArray,
                    eccentricity: FloatArray, obliquity: FloatArray,
                    semi_major_axis: FloatArray, radius: float,
                    use_obliquity: bool = True, eccentricity_truncation_lvl: int = 2,
                    max_order_l: int = 2) -> Tuple[Dict[FreqSig, FloatArray], Dict[FreqSig, Dict[int, DissipTermsMix]]]:
    """ Calculate tidal dissipation terms for each tidal harmonic order l

    Returns the tidal heating and tidal potential derivative components for each unique frequency. They still must be
        multiplied by the Love number (the Love's sign will be taken into account here and should not be
        multiplied later), added together, and then scaled by the tidal susceptibility.

    This function wraps tides.calculate_terms_at_orderl which is njit'd. This function can not be njit'd with the
        current version of numba (0.46) due to the different eccentricity/inclination functions. Even though these
        functions are njit'd, there is not an njit-safe way to access them within an njit function (without compiling
        them all).

    Parameters
    ----------
    orbital_frequency : FloatArray
        Planet's orbital frequency around tidal host [rad s-1]
    spin_frequency : FloatArray
        Planet's rotation frequency [rad s-1]
    eccentricity : FloatArray
        Planet's orbital eccentricity around tidal host
    obliquity : FloatArray
        Planet's axial tilt relative to orbital plane around tidal host [radians]
    semi_major_axis : FloatArray
        Orbital separation [m]
    radius : float
        Planet's radius (used to calculate R/a) [m]
    use_obliquity : bool = True
        Determine if a non-zero obliquity is allowed.
        If you are sure that obliquity should always be equal to zero, then set to this False for faster computation.
    eccentricity_truncation_lvl : int = 2
        Eccentricity functions must be truncated to a power of e. The higher the truncation, the greater the accuracy
        of the result, especial when dealing with high eccentricities.
        Truncation level 2 is recommended when e <~ 0.05. See Renaud et al. 2020 for more information.
    max_order_l : int = 2
        Maximum harmonic number to include

    Returns
    -------
    unique_freqs : Dict[FreqSig, FloatArray]
        Each unique frequency stored as a signature (orbital motion and spin-rate coeffs), and the calculated frequency
            (combo of orbital motion and spin-rate) [rad s-1]
    results_by_unique_frequency : Dict[FreqSig, Dict[int, DissipTermsMix]]
        Results for tidal heating, dU/dM, dU/dw, dU/dO are stored in a tuple for each tidal harmonic l and
            unique frequency.

    See Also
    --------
    TidalPy.tides.dissipation.calculate_terms_at_orderl
    """

    # Get eccentricity functions at desired truncation level (multiple functions for each order-l).
    eccentricity_funcs_by_orderl = eccentricity_truncations[eccentricity_truncation_lvl]

    # Get inclination (obliquity) functions based on *if* obliquity is always equal to zero or not (multiple functions
    #    for each order-l).
    inclination_funcs_by_orderl = inclination_functions[use_obliquity]

    # Store results by each unique frequency so that outer-scope funcs can parse through the unique frequencies and
    #    calculate the complex compliances only once.
    results_by_unique_frequency = dict()
    unique_freqs = dict()

    # This for loop is generally not going to be large (for example, it will only have one loop for max_l == 2)
    for order_l in range(2, max_order_l + 1):
        # Pull out the l-specific eccen/inclin funcs and calculate them.
        eccentricity_results_at_orderl = eccentricity_funcs_by_orderl[order_l](eccentricity)
        obliquity_results_at_orderl = inclination_funcs_by_orderl[order_l](obliquity)

        # Sort universal coefficients (multipliers to all tidal terms) by order l
        universal_coeffs_by_m = get_universal_coeffs(order_l)

        unique_freqs_at_this_l, results_at_this_l_by_frequency = \
            calculate_terms_at_orderl(orbital_frequency, spin_frequency, semi_major_axis, radius,
                                      eccentricity_results_at_orderl, obliquity_results_at_orderl,
                                      universal_coeffs_by_m, order_l)

        # Combine all order-l unique frequencies into a global unique freq list.
        for (freq_sig, freq), dissipation_terms in zip(unique_freqs_at_this_l, results_at_this_l_by_frequency):
            if freq_sig not in unique_freqs:
                unique_freqs[freq_sig] = freq
                results_by_unique_frequency[freq_sig] = {order_l: dissipation_terms}
            else:
                results_by_unique_frequency[freq_sig][order_l] = dissipation_terms

    return unique_freqs, results_by_unique_frequency

@njit
def calculate_terms_at_orderl(orbital_frequency: FloatArray, spin_frequency: FloatArray,
                              semi_major_axis: FloatArray, radius: float,
                              eccentricity_results: Dict[int, Dict[int, FloatArray]],
                              inclination_results: Dict[Tuple[int, int], FloatArray],
                              universal_coeff_by_m: Dict[int, float],
                              order_l: int) -> Tuple[List[Tuple[FreqSig, FloatArray]], List[DissipTermsMix]]:
    """ Calculate tidal dissipation terms for a given tidal harmonic order-l.

    Requires the user to provide the correct eccentricity and inclination functions. This is generally done by the
    tides.calculate_terms wrapper.

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Planet's orbital frequency around tidal host [rad s-1]
    spin_frequency : np.ndarray
        Planet's rotation frequency [rad s-1]
    semi_major_axis : np.ndarray
        Orbital separation [m]
    radius : float
        Planet's radius (used to calculate R/a) [m]
    eccentricity_results : Dict[int, Dict[int, FloatArray]]
        Eccentricity results pre-calculated using the tides.eccentricityFuncs functions.
        Stored in a dict(p: dict(q: result))
    inclination_results : Dict[Tuple[int, int], FloatArray]
        Inclination/Obliquity results pre-calculated using the tides.inclinationFuncs functions.
        Stored in a dict(tuple(m, p): result)
    universal_coeff_by_m : Dict[int, float]
        Universal coefficients for this order-l.
        Stored in a dict(m: result)
    order_l : int
        Tidal harmonic order-l

    Returns
    -------
    unique_frequencies : List[Tuple[FreqSig, FloatArray]]
        Unique frequency signature and calculate frequency [rad s-1]
        Frequencies are integer combinations of orbital motion and spin-rate.
    results_by_frequency : List[DissipTermsMix]
        Tidal dissipation terms stored for each unique frequency.
        At each frequency dissipation terms are stored as:
            heating_term, dUdM_term, dUdw_term, dUdO_term
            These still need to be multiplied by the -Im[k_l] and tidal susceptibility.

    See Also
    --------
    TidalPy.tides.dissipation.calculate_terms
    """

    # Storage for results by unique frequency
    results_by_frequency = list()
    unique_frequencies = list()
    encountered_sigs = list()

    # Order l controls the distance multiplier. The minus 4 comes from the tidal susceptibility already carrying
    #    (R / a)^5 so (R / a)^(2l + 1) --> (R / a)^(2l - 4)
    distance_multiplier = (radius / semi_major_axis)**(2 * order_l - 4)

    for (m, p), obliquity_terms in inclination_results.items():
        # Pull out the m and p integers from the non-zero obliquity terms

        # Pull out universal coefficient
        #    The Tidal Susceptibility carries a factor of (3 / 2) already. So we need to divide the uni_coeff
        #    by that much to ensure it is not double counted.
        uni_coeff = distance_multiplier * universal_coeff_by_m[m] / 1.5

        for q, eccentricity_terms in eccentricity_results[p].items():
            # Pull out q integer from the non-zero eccentricity terms

            # Multiply eccentricity terms by the obliquity terms
            eccen_obliq = obliquity_terms * eccentricity_terms

            # Calculate tidal mode, frequency, and sign
            n_coeff = (order_l - 2 * p + q)

            # Skip terms that we know are zero (0 * n - 0 * m will always be 0 freq and -imk(0) == 0)
            if n_coeff == 0 and m == 0:
                continue

            n_coeff_abs = abs(n_coeff)
            mode = n_coeff * orbital_frequency - m * spin_frequency
            mode_sign = np.sign(mode)
            mode_frequency = np.abs(mode)

            # Determine the frequency signature used to store unique frequencies
            #    FIXME: Note: numba does not support f-strings currently
            if n_coeff == 0:
                freq_sig = (0, m)
            elif m == 0:
                freq_sig = (n_coeff_abs, 0)
            elif n_coeff < 0:
                freq_sig = (n_coeff_abs, m)
            else:
                freq_sig = (n_coeff, -m)

            # Calculate coefficients for heating and potential derivatives
            uni_multiplier = uni_coeff * eccen_obliq
            heating_term = uni_multiplier * mode_frequency
            dUdM_term = uni_multiplier * n_coeff * mode_sign
            dUdw_term = uni_multiplier * (order_l - 2. * p) * mode_sign
            dUdO_term = uni_multiplier * m * mode_sign

            # The tidal heating and potential derivatives should also be multiplied by the love number calculated
            #    for that mode. But, the Love number only cares about the complex compliance which in turn only requires
            #    frequencies (does not care about modes or order-l). Therefore, we can store everything in terms of
            #    *unique* frequencies.

            # TODO: for now njit really does not like dicts. Even typed dict did not seem to work here
            #    (as of numba v.46). So, instead, we will perform a loop to look at frequency signatures and only record
            #    modes when a unique signature is found.

            i = 0
            found = False
            found_at = 0
            for freq_sig_test in encountered_sigs:
                if freq_sig == freq_sig_test:
                    # Once (if) a duplicate signature is found break out of this loop.
                    found = True
                    found_at = i
                    break
                i += 1

            if found:
                # For previously found signatures we need to add the tidal dissipation terms to the ones already
                #    recorded.
                heating_term_old, dUdM_term_old, dUdw_term_old, dUdO_term_old = results_by_frequency[found_at]

                heating_combo = heating_term + heating_term_old
                dUdM_combo = dUdM_term + dUdM_term_old
                dUdw_combo = dUdw_term + dUdw_term_old
                dUdO_combo = dUdO_term + dUdO_term_old

                # Replace the output with the newly summed value. This should work for floats and arrays.
                results_by_frequency[found_at] = (heating_combo, dUdM_combo, dUdw_combo, dUdO_combo)
            else:
                # If the frequency signature is not found, make a new record for this frequency.
                results_by_frequency.append((heating_term, dUdM_term, dUdw_term, dUdO_term))
                unique_frequencies.append((freq_sig, mode_frequency))
                encountered_sigs.append(freq_sig)

    return unique_frequencies, results_by_frequency

@njit
def mode_collapse(gravity: float, radius: float, density: float, shear_modulus: FloatArray,
                  complex_compliance_by_frequency: Dict[str, ComplexArray],
                  tidal_terms_by_frequency: Dict[FreqSig, Dict[int, DissipTermsMix]],
                  tidal_susceptibility: FloatArray, tidal_host_mass: float, tidal_scale: float,
                  max_tidal_l: int) -> Tuple[FloatArray, FloatArray, FloatArray, FloatArray, ComplexArray, FloatArray]:
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
    complex_compliance_by_frequency : Dict[str, ComplexArray]
        The complex compliance for the layer or planet calculated at each unique tidal frequency [Pa-1]
    tidal_terms_by_frequency : Dict[FreqSig, Dict[int, DissipTerms]]
        Each dissipation term: E^dot, dUdM, dUdw, dUdO; calculated for each unique tidal frequency and order-l
    tidal_susceptibility : FloatArray
        Tidal susceptibility, defined as (3/2) G M_host^2 R_target^5 / a^6 [N m]
    tidal_host_mass : float
        Mass of the tidal host [kg]
        Needed to offset the tidal susceptibility for the potential derivatives.
    tidal_scale : float
        Multiplier that ranges from 0 to 1 that scales the imaginary portion of the complex love number to simulate
            only a portion of a planet's volume is contributing to the tidal dissipation.
    max_tidal_l : int
        Tidal harmonic order

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

    tidal_heating_terms = list()
    dUdM_terms = list()
    dUdw_terms = list()
    dUdO_terms = list()
    love_number_terms = list()
    negative_imk_terms = list()

    for tidal_order_l in range(2, max_tidal_l + 1):

        # Calculate the effective rigidity which does not use the complex compliance.
        # TODO: Should this use the static shear modulus or the effective (post-melting/heating) shear?
        #    Right now it is using the effective rigidity. If it only needs the static then this could be done outside
        #    this function.
        effective_rigidity = effective_rigidity_general(shear_modulus, gravity, radius, density,
                                                        order_l=tidal_order_l)

        # Pull out the already computed complex compliances for each frequency
        for unique_freq_signature, complex_compliance in complex_compliance_by_frequency.items():

            # Calculate the complex Love number, which also changes functional form for each order_l
            complex_love = complex_love_general(complex_compliance, shear_modulus, effective_rigidity,
                                                order_l=tidal_order_l)

            # Scale the Love number by the layer's contribution
            # TODO: Should the tidal scale affect the Re(k) as well as the Im(k)?
            complex_love = np.real(complex_love) + (np.imag(complex_love) * tidal_scale)
            neg_imk = -np.imag(complex_love)
            neg_imk_scaled = neg_imk * tidal_susceptibility

            # The tidal potential carries one fewer dependence on the tidal host mass that is built into the
            #    tidal susceptibility. Divide that out now.
            neg_imk_scaled_potential = neg_imk_scaled / tidal_host_mass

            # Pull out the tidal terms pre-calculated for this unique frequency. See Tides.orbital_change
            heating_term, dUdM_term, dUdw_term, dUdO_term = \
                tidal_terms_by_frequency[unique_freq_signature][tidal_order_l]

            # Store results
            tidal_heating_terms.append(heating_term * neg_imk_scaled)
            dUdM_terms.append(dUdM_term * neg_imk_scaled_potential)
            dUdw_terms.append(dUdw_term * neg_imk_scaled_potential)
            dUdO_terms.append(dUdO_term * neg_imk_scaled_potential)
            love_number_terms.append(complex_love)
            negative_imk_terms.append(neg_imk_scaled)

    # Collapse Modes
    tidal_heating = sum(tidal_heating_terms)
    dUdM = sum(dUdM_terms)
    dUdw = sum(dUdw_terms)
    dUdO = sum(dUdO_terms)
    love_number = sum(love_number_terms)
    negative_imk = sum(negative_imk_terms)

    return tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk
