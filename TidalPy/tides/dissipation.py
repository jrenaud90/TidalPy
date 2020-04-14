from typing import Tuple, Dict

import numpy as np
from numba import njit
from scipy.constants import G

from ..types import FloatArray
from .universal_coeffs import get_universal_coeffs
from .inclinationFuncs import inclination_functions
from .eccentricityFuncs import eccentricity_truncations


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

@njit
def mode_collapse( orbital_frequency: FloatArray, spin_frequency: FloatArray,
                  eccentricity: FloatArray, obliquity: FloatArray,
                  use_obliquity: bool = True, eccentricity_truncation_lvl: int = 2, max_order_l: int = 2) -> \
        Dict[str, Tuple[FloatArray, Dict[int, Tuple[FloatArray, FloatArray, FloatArray, FloatArray]]]]:
    """ Collapses tidal modes down to their unique frequencies for tidal heating and the tidal potential derivatives.

    Returns the tidal heating and tidal potential derivative components for each unique frequency. They still must be
        multiplied by the Love number (the Love's sign will be taken into account here and should not be
        multiplied later), added together, and then scaled by the tidal susceptibility.

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Planet's orbital frequency around tidal host [rad s-1]
    spin_frequency : np.ndarray
        Planet's rotation frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's orbital eccentricity around tidal host
    obliquity : np.ndarray
        Planet's axial tilt relative to orbital plane around tidal host [radians]
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
    results_by_uniquefreq : Dict[str, Tuple[FloatArray, Dict[int, Tuple[FloatArray, FloatArray, FloatArray, FloatArray]]]]
        Results for tidal heating, dU/dM, dU/dw, dU/dO are stored in a tuple for each tidal harmonic l and
        unique frequency.
    """

    # Determine if obliquity is used
    obliquity_func_set = inclination_functions[use_obliquity]
    eccentricity_func_set = eccentricity_truncations[eccentricity_truncation_lvl]

    # Storage for results by unique frequency
    results_by_uniquefreq = dict()

    # This for loop is generally not going to be large (often will only loop once, for max_l == 2)
    for order_l in range(2, max_order_l + 1):

        # Find the obliquity and eccentricity functions for this harmonic
        eccentricity_at_orderl = eccentricity_func_set[order_l](eccentricity)
        obliquity_at_orderl = obliquity_func_set[order_l](obliquity)

        # Sort universal coefficients (multipliers to all tidal terms) by order l
        universal_coeffs_bym = get_universal_coeffs(order_l)

        for (p, m), obliquity_terms in obliquity_at_orderl.items():
            # Pull out the p and m integers from the non-zero obliquity terms

            # Pull out universal coefficient
            #    The Tidal Susceptibility carries a factor of (3 / 2) already. So we need to divide the uni_coeff
            #    by that much to ensure it is not double counted.
            uni_coeff = universal_coeffs_bym[m] / 1.5

            for q, eccentricity_terms in eccentricity_at_orderl[p].items():
                # Pull out q integer from the non-zero eccentricity terms

                # Multiply eccentricity terms by the obliquity terms
                eccen_obliq = obliquity_terms * eccentricity_terms

                # Calculate tidal mode, frequency, and sign
                n_coeff = (order_l - 2 * p + q)
                mode = n_coeff * orbital_frequency - m * spin_frequency
                mode_sign = np.sign(mode)
                mode_freq = np.abs(mode)

                # Determine the frequency signature used to store unique frequencies
                if n_coeff == 0:
                    freq_sig = f'{m}|S|'
                elif n_coeff < 0:
                    freq_sig = f'|{n_coeff}n + {m}S|'
                else:
                    freq_sig = f'|{n_coeff}n - {m}S|'

                # Calculate coefficients for heating and potential derivatives
                uni_multiplier = uni_coeff * eccen_obliq
                heating_term = uni_multiplier * mode_freq
                dUdM_term = uni_multiplier * n_coeff * mode_sign
                dUdw_term = uni_multiplier * (order_l - 2. * p) * mode_sign
                dUdO_term = uni_multiplier * m * mode_sign

                # The tidal heating and potential derivatives should also be multiplied by the love number calculated
                #    for that mode. But, love number only cares about the frequencies and some of them may be repeated

                if freq_sig not in results_by_uniquefreq:
                    # New frequency - create dictionary to store different order-l results.
                    results_by_uniquefreq[freq_sig] = (mode_freq, dict())
                    results_by_uniquefreq[freq_sig][1][order_l] = (heating_term, dUdM_term, dUdw_term, dUdO_term)
                else:
                    if order_l in results_by_uniquefreq[freq_sig][1]:
                        results_by_uniquefreq[freq_sig][1][order_l][0] += heating_term
                        results_by_uniquefreq[freq_sig][1][order_l][1] += dUdM_term
                        results_by_uniquefreq[freq_sig][1][order_l][2] += dUdw_term
                        results_by_uniquefreq[freq_sig][1][order_l][3] += dUdO_term
                    else:
                        results_by_uniquefreq[freq_sig][1][order_l] = (heating_term, dUdM_term, dUdw_term, dUdO_term)

    return results_by_uniquefreq


