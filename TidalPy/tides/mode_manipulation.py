from typing import List, Tuple, Dict, Union

import numpy as np

from .universal_coeffs import get_universal_coeffs
from .inclinationFuncs import inclination_functions
from .eccentricityFuncs import eccentricity_truncations
from .love1d import effective_rigidity_general, complex_love_general
from ..utilities.performance.numba import njit
from ..utilities.types import FloatArray, ComplexArray, NoneType

FreqSig = Tuple[int, int]
DissipTermsMix = Tuple[FloatArray, FloatArray, FloatArray, FloatArray]

def build_mode_manipulators(max_order_l: int = 2, eccentricity_truncation_lvl: int = 8, use_obliquity: bool = True):

    # Pull out eccentricity and obliquity (inclination) functions for each order-l.
    #    Note for all lists and tuples that are of order-l, assume that they are off by 2:
    #    l=2 -> index=0, l=3 -> index=1
    obliquity_func_by_order_l_temp = inclination_functions[use_obliquity]
    eccentricity_func_by_order_l_temp = eccentricity_truncations[eccentricity_truncation_lvl]

    # Store functions in list (that will then be tupled) this will interface in a numba-safe way.
    obliquity_func_by_order_l_minus2 = list()
    eccentricity_func_by_order_l_minus2 = list()
    for order_l_ in range(2, max_order_l + 1):
        obliquity_func_by_order_l_minus2.append(obliquity_func_by_order_l_temp[order_l_])
        eccentricity_func_by_order_l_minus2.append(eccentricity_func_by_order_l_temp[order_l_])
    obliquity_func_by_order_l_minus2 = tuple(obliquity_func_by_order_l_minus2)
    eccentricity_func_by_order_l_minus2 = tuple(eccentricity_func_by_order_l_minus2)

    # Since this name space is used for the following functions, let's make sure it is as clean as possible
    del eccentricity_func_by_order_l_temp, obliquity_func_by_order_l_temp, eccentricity_truncation_lvl,\
        use_obliquity, order_l_

    @njit
    def calculate_terms(orbital_frequency: FloatArray, spin_frequency: FloatArray,
                        eccentricity: FloatArray, obliquity: FloatArray,
                        semi_major_axis: FloatArray, radius: float) \
            -> Tuple[Dict[FreqSig, FloatArray], Dict[FreqSig, Dict[int, DissipTermsMix]]]:
        """ Calculate tidal dissipation terms and frequencies based on the current orbital and spin state.

        Requires the user to provide the correct eccentricity and inclination functions. This is generally done by the
        tides.calculate_terms wrapper.

        Parameters
        ----------
        orbital_frequency : FloatArray
            Orbit's orbital frequency around tidal host [rad s-1]
        spin_frequency : FloatArray
            World's rotation frequency [rad s-1]
        eccentricity : FloatArray
            Orbit's eccentricity (must be closed orbits with 0 <= e < 1)
        obliquity : FloatArray
            World's axial tilt relative to the orbital plane [radians]
        semi_major_axis : FloatArray
            Orbital separation [m]
        radius : float
            Planet's radius (used to calculate R/a) [m]


            # TODO
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

        # Storage for results by unique frequency signature. Must provide fake values in order to initialize
        #    numba.
        fake_result = orbital_frequency * spin_frequency * eccentricity * obliquity
        results_by_frequency = {
            (-100, -100): {
                0: fake_result
            }
        }
        unique_frequencies = {
            (-100, -100): orbital_frequency - spin_frequency
        }

        for order_l in range(2, max_order_l + 1):

            # Get universal coefficients for this order l
            universal_coeff_by_m = get_universal_coeffs(order_l)

            # Pull out eccentricity and obliquity function (from global variable stroage), and calculate
            eccentricity_func = eccentricity_func_by_order_l_minus2[order_l - 2]
            obliquity_func = obliquity_func_by_order_l_minus2[order_l - 2]
            eccentricity_results = eccentricity_func(eccentricity)
            obliquity_results = obliquity_func(obliquity)

            # Order l controls the distance multiplier. The minus 4 comes from the tidal susceptibility already carrying
            #    (R / a)^5 so (R / a)^(2l + 1) --> (R / a)^(2l - 4)
            distance_multiplier = (radius / semi_major_axis)**(2 * order_l - 4)

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
                    if freq_sig not in unique_frequencies:
                        # New unique frequency
                        unique_frequencies[freq_sig] = mode_frequency
                        results_by_frequency[freq_sig] = {order_l: (heating_term, dUdM_term, dUdw_term, dUdO_term)}
                    else:
                        if order_l in results_by_frequency[freq_sig]:
                            # Previous results found at this frequency and order_l, add new results to previous one
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


    def collapse_modes(gravity: float, radius: float, density: float, shear_modulus: Union[NoneType, FloatArray],
                       complex_compliance_by_frequency: Dict[FreqSig, ComplexArray],
                       tidal_terms_by_frequency: Dict[FreqSig, Dict[int, DissipTermsMix]],
                       tidal_susceptibility: FloatArray, tidal_host_mass: float, tidal_scale: float,
                       cpl_ctl_method: bool = False) -> \
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
        complex_compliance_by_frequency : Dict[FreqSig, ComplexArray]
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

        tidal_heating_terms = list()
        dUdM_terms = list()
        dUdw_terms = list()
        dUdO_terms = list()
        love_number_terms = list()
        negative_imk_terms = list()

        for tidal_order_l in range(2, max_order_l + 1):

            # Calculate the effective rigidity which does not use the complex compliance.
            # TODO: Should this use the static shear modulus or the effective (post-melting/heating) shear?
            #    Right now it is using the effective rigidity. If it only needs the static then this could be done outside
            #    this function.
            if cpl_ctl_method:
                effective_rigidity = shear_modulus
            else:
                if shear_modulus is None:
                    raise Exception('Shear modulus is required for non-CTL/CPL models. Set it to a fake float e.g., 1.')

                effective_rigidity = effective_rigidity_general(shear_modulus, gravity, radius, density,
                                                                order_l=tidal_order_l)

            # Pull out the already computed complex compliances for each frequency
            for unique_freq_signature, complex_compliance in complex_compliance_by_frequency.items():

                if cpl_ctl_method:
                    # In the CTL/CPL method, the complex love number is passed in as the complex compliance
                    complex_love = complex_compliance
                else:
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

                # Pull out the tidal terms pre-calculated for this unique frequency. See Tides.update_orbit_spin
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


    return calculate_terms, collapse_modes