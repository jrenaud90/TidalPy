import numpy as np
from numba import njit

from .love_1d import effective_rigidity_general, complex_love_general
from .universal_coeffs import get_universal_coeffs
from .inclinationFuncs import inclination_functions
from .eccentricityFuncs import eccentricity_truncations
from ..rheology.compliance_models import fixed_q


# @njit
def kaula_collapse(spin_frequency, orbital_frequency, semi_major_axis,
                   eccentricity_results_byorderl, inclination_results_byorderl,
                   complex_compliance_func, complex_compliance_input,
                   shear_modulus, viscosity, planet_radius, planet_gravity, planet_density,
                   max_order_l):
    """
    """

    compliance = 1. / shear_modulus
    cached_complex_compliances = {}
    cached_freqs = {}

    tidal_modes = []
    tidal_heating_bymode = []
    dUdM_bymode = []
    dUdw_bymode = []
    dUdO_bymode = []

    sigs = []

    for order_l in range(2, max_order_l+1):

        effective_rigidity = effective_rigidity_general(shear_modulus, planet_gravity, planet_radius, planet_density,
                                                        order_l=order_l)

        inclination_results = inclination_results_byorderl[order_l-2]
        eccentricity_results = eccentricity_results_byorderl[order_l-2]
        universal_coeffs = get_universal_coeffs(order_l-2)

        # Order l acts as an exponent to the radius / semi-major axis. The tidal susceptibility already considers l=2
        #    It contains R^5 / a^6 already so that is why we subtract a 4 off below (2l + 1) @ l=2.
        if order_l > 2:
            distance_scale = (planet_radius / semi_major_axis)**(2 * order_l - 4)
        else:
            distance_scale = np.ones_like(semi_major_axis)

        for (_, p, q), eccentricity_result in eccentricity_results.items():

            # Pull out specific inclination result
            inclination_subresults = inclination_results[p]

            for m in range(0, order_l + 1):

                # Driving term: eccentricity function squared times inclination function squared
                F2G2 = eccentricity_result * inclination_subresults[m]

                # Find universal coefficient.
                #    The Tidal Susceptibility carries a factor of (3 / 2) already. So we need to divide the uni_coeff
                #    by that much to ensure it is not double counted.
                uni_coeff = universal_coeffs[m] / 1.5

                # Collapse all multipliers
                multiplier = distance_scale * uni_coeff * F2G2

                # Find tidal mode and frequency
                orbital_coeff = (order_l - 2 * p + q)
                spin_coeff = -m
                mode = orbital_coeff * orbital_frequency + spin_coeff * spin_frequency
                freq = np.abs(mode)
                sgn = np.sign(mode)

                # Complex compliances are only a function of frequency (at least at this stage)
                #   Since they can be expensive to calculate, and since the number of unique frequencies is < #modes
                #   We will use a frequency signature to store cached complex compliance results
                #   This is especially important for max_order_l > 2 as many duplicate freqs will be hit.
                mode_signature = (orbital_coeff, spin_coeff)
                if mode_signature in cached_complex_compliances:
                    complex_compliance = cached_complex_compliances[mode_signature]
                else:
                    complex_compliance = complex_compliance_func(compliance, viscosity, freq, *complex_compliance_input)
                    cached_complex_compliances[mode_signature] = complex_compliance
                    cached_freqs[mode_signature] = freq

                if complex_compliance_func is fixed_q:
                    Q, _, k_2 = complex_compliance_input
                    cmplx_love_number = k_2 * (1.0 - 1.0j / Q) * np.ones_like(freq * shear_modulus)
                else:
                    cmplx_love_number = complex_love_general(complex_compliance, shear_modulus, effective_rigidity, order_l=order_l)

                neg_imk = -np.imag(cmplx_love_number)
                neg_imk_sgn = sgn * neg_imk

                # Store Results
                sigs.append(mode_signature)

                tidal_modes.append(mode)
                tidal_heating_bymode.append(freq * neg_imk * multiplier)
                dUdM_bymode.append(orbital_coeff * neg_imk_sgn * multiplier)
                dUdw_bymode.append((order_l - 2*p) * neg_imk_sgn * multiplier)
                dUdO_bymode.append(m * neg_imk_sgn * multiplier)


    if complex_compliance_func is not fixed_q:
        print(complex_compliance_func)
        from pprint import pprint

        heatings = {}
        for sig, heating in zip(sigs, tidal_heating_bymode):
            if sig in heatings:
                heatings[sig] += heating
            else:
                heatings[sig] = heating
        pprint(heatings)
        pprint([(mode, tidal_heating_) for mode, tidal_heating_ in zip(sigs, dUdM_bymode)])
        pprint([(mode, tidal_heating_) for mode, tidal_heating_ in zip(sigs, dUdw_bymode)])
        pprint([(mode, tidal_heating_) for mode, tidal_heating_ in zip(sigs, dUdO_bymode)])
    return tidal_modes, tidal_heating_bymode, dUdM_bymode, dUdw_bymode, dUdO_bymode


def calculate(spin_frequency, orbital_frequency, semi_major_axis, eccentricity, inclination,
              complex_compliance_func, complex_compliance_input,
              shear_modulus, viscosity, planet_radius, planet_gravity, planet_density,
              eccentricity_truncation: int = 6, use_inclination: bool = True, max_order_l: int = 2):
    """ Calculate tidal potential derivatives and tidal heating """

    # Check if inclination was provided as an input
    if inclination is None:
        use_inclination = False

    # Get the truncation level for eccentricity
    eccentricity_func_byorder = eccentricity_truncations[eccentricity_truncation]

    # See if inclination is being used. If so then grab the appropriate function
    inclination_func_byorder = inclination_functions[use_inclination]

    eccentricity_results_byorderl = list()
    inclination_results_byorderl = list()
    for order_l in range(2, max_order_l+1):

        # Get eccentricity function calculators
        eccentricity_functions_sqrd = eccentricity_func_byorder[order_l-2](eccentricity)
        eccentricity_results_byorderl.append(eccentricity_functions_sqrd)

        # Get inclination function calculators
        inclination_functions_sqrt = inclination_func_byorder[order_l-2](inclination)
        inclination_results_byorderl.append(inclination_functions_sqrt)

    # Calculate heating and potential derivative coefficients
    tidal_modes, tidal_heating_bymode, dUdM_bymode, dUdw_bymode, dUdO_bymode = \
        kaula_collapse(spin_frequency, orbital_frequency, semi_major_axis,
                       eccentricity_results_byorderl, inclination_results_byorderl,
                       complex_compliance_func, complex_compliance_input,
                       shear_modulus, viscosity, planet_radius, planet_gravity, planet_density,
                       max_order_l)

    # Collapse modes down to single values
    tidal_heating = sum(tidal_heating_bymode)
    dUdM = sum(dUdM_bymode)
    dUdw = sum(dUdw_bymode)
    dUdO = sum(dUdO_bymode)

    return tidal_heating, dUdM, dUdw, dUdO



