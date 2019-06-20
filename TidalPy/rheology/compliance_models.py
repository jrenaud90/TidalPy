""" Rheology Compliance Functions

How To Implement a New Rheology:
    * Make a new function here with a function doc-string that follows the format of the pre-built rheologies.
        - Items under the 'Parameters' header will be used to help construct the rheology. You need to supply them.
        - All rheologies require three arguments: compliance (non-complex), viscosity, and frequency.
            These should be the first arguments and in this order.
        - Additional arguments can be supplied by adding them after these required arguments but they must be included in the function's
            doc-string after the 'other args' key (in the same order).
    * Any new arguments that are simple constants can then be added to the rheology section of a planet's layer configuration.
        - Non-constants will require additional coding. As of right now there is not a clean way to add such functionality. But, you can
            look at the zeta-frequency dependency implementation.
"""

import numpy as np

from ..performance import njit, find_factorial

@njit
def off(compliance, viscosity, frequency):
    """ No Rheology

    other args: None
    """

    return compliance + 0.0j


@njit
def fixed_q(compliance, viscosity, frequency, quality_factor, planet_beta):
    """ Fixed-Q Rheology

        --- Parameters ---
        nice name:  Fixed-Q
        line style: -.
        color:      k
        other args: quality_factor, planet_beta

    """

    real_j = -19. / (2. * planet_beta)
    imag_j = -quality_factor * (19. / 4.) * (2. * planet_beta * compliance + 19.) / (compliance * planet_beta**2)
    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def maxwell(compliance, viscosity, frequency):
    """ Maxwell Rheology

        --- Parameters ---
        nice name:  Maxwell
        line style: -
        color:      k
        other args: None

    """

    real_j = compliance
    imag_j = -1.0 / (viscosity * frequency)

    complex_compliance = real_j + 1.0j*imag_j

    return complex_compliance


@njit
def voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset):
    """ Voigt-Kelvin Rheology

        --- Parameters ---
        nice name:  Voigt-Kelvin
        line style: -
        color:      g
        other args: voigt_compliance_offset, voigt_viscosity_offset

    """

    voigt_comp = voigt_compliance_offset * compliance
    voigt_visc = voigt_viscosity_offset * viscosity

    denominator = (voigt_comp * voigt_visc * frequency)**2 + 1.
    real_j = voigt_comp / denominator
    imag_j = -voigt_comp**2 * voigt_visc * frequency / denominator
    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance


@njit
def burgers(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset):
    """ Burgers Rheology

        --- Parameters ---
        nice name:  Burgers
        line style: -
        color:      r
        other args: voigt_compliance_offset, voigt_viscosity_offset

    """

    voigt_complex_comp = voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset)
    maxwell_complex_comp = maxwell(compliance, viscosity, frequency)

    return voigt_complex_comp + maxwell_complex_comp


@njit
def andrade(compliance, viscosity, frequency, alpha, zeta):
    """ Andrade Rheology

        --- Parameters ---
        nice name:  Andrade
        line style: -
        color:      b
        other args: alpha, zeta

    """
    maxwell_complex_comp = maxwell(compliance, viscosity, frequency)

    const_term = compliance * (compliance * viscosity * frequency * zeta)**(-alpha) * \
                               find_factorial(alpha)
    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term
    andrade_complex_comp = real_j + 1.0j*imag_j

    return maxwell_complex_comp + andrade_complex_comp


@njit
def andrade_freq(compliance, viscosity, frequency, alpha, zeta, andrade_freq_params, andrade_freq_func):
    """ Andrade Rheology with frequency-dependent Zeta

        --- Parameters ---
        nice name:  Andrade
        line style: -
        color:      b
        other args: alpha, zeta, andrade_freq_params, andrade_freq_func

    """
    # The andrade_freq_func is defined by one of the models in andrade_frequency.py
    zeta = andrade_freq_func(zeta, frequency, *andrade_freq_params)

    maxwell_complex_comp = maxwell(compliance, viscosity, frequency)

    const_term = compliance * (compliance * viscosity * frequency * zeta)**(-alpha) * \
                               find_factorial(alpha)
    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term
    andrade_complex_comp = real_j + 1.0j * imag_j

    return maxwell_complex_comp + andrade_complex_comp


@njit
def sundberg(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta):
    """ Sundberg-Cooper Rheology

        --- Parameters ---
        nice name:  Sundberg
        line style: -
        color:      m
        other args: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta

    """
    voigt_complex_comp = voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset)
    andrade_complex_comp = andrade(compliance, viscosity, frequency, alpha, zeta)

    return voigt_complex_comp + andrade_complex_comp


@njit
def sundberg_freq(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta,
                  andrade_freq_params, andrade_freq_func):
    """ Sundberg-Cooper Rheology

        --- Parameters ---
        nice name:  Sundberg
        line style: -
        color:      m
        other args: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta, andrade_freq_params, andrade_freq_func

    """
    voigt_complex_comp = voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset)
    andrade_complex_comp = andrade_freq(compliance, viscosity, frequency, alpha, zeta,
                                        andrade_freq_params, andrade_freq_func)

    return voigt_complex_comp + andrade_complex_comp