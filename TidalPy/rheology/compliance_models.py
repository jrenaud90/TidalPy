""" Complex Compliance Functions

    How To Implement a New Rheological Model:

        * Make a new function here with a function doc-string that follows the format of the pre-built rheologies.
            - Items under the 'Parameters' header will be used to help construct the rheology. You need to supply them.
            - All rheologies require three arguments: compliance (non-complex), viscosity, and frequency.
                These should be the first arguments and in this order.
            - Additional arguments can be supplied by adding them after these required arguments but they must be
                included in the function's doc-string after the 'other args' key (in the same order).

        * Any new arguments that are simple constants can then be added to the rheology section of a planet's
            layer configuration.
            - Non-constants will require additional coding. As of right now there is not a clean way to add such
                functionality. But, you can look at the zeta-frequency dependency implementation.

"""

import numpy as np

from ..performance import find_factorial, tpy_vectorize
from ..types import float_eps


@tpy_vectorize(['complex128(float64, float64, float64)'])
def off(compliance, viscosity, frequency):
    """ No Tides

    other args: None
    """

    real_j = compliance
    imag_j = 0.

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance


@tpy_vectorize(['complex128(float64, float64, float64, float64, float64)'])
def fixed_q(compliance, viscosity, frequency, quality_factor, planet_beta):
    """ Fixed-Q Tides

    !TPY_args const: quality_factor, planet_beta

    """

    real_j = -19. / (2. * planet_beta)
    if np.abs(frequency) <= float_eps:
        imag_j = 0.
    else:
        imag_j = -quality_factor * (19. / 4.) * (2. * planet_beta * compliance + 19.) / (compliance * planet_beta**2) \
            * (1. + 0. * frequency)

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance


@tpy_vectorize(['complex128(float64, float64, float64)'])
def maxwell(compliance, viscosity, frequency):
    """ Maxwell Tides

    !TPY_args const: None

    """

    real_j = compliance
    if np.abs(frequency) <= float_eps:
        imag_j = 0.
    else:
        imag_j = -1. / (viscosity * frequency)

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance


@tpy_vectorize(['complex128(float64, float64, float64, float64, float64)'])
def voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset):
    """ Voigt-Kelvin Tides

    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    """

    voigt_comp = voigt_compliance_offset * compliance
    voigt_visc = voigt_viscosity_offset * viscosity

    if np.abs(frequency) <= float_eps:
        real_j = 0.
        imag_j = 0.
    else:
        denominator = (voigt_comp * voigt_visc * frequency)**2 + 1.
        real_j = voigt_comp / denominator
        imag_j = -voigt_comp**2 * voigt_visc * frequency / denominator

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance


@tpy_vectorize(['complex128(float64, float64, float64, float64, float64)'])
def burgers(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset):
    """ Burgers Tides

    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    """

    voigt_complex_comp = voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset)
    maxwell_complex_comp = maxwell(compliance, viscosity, frequency)

    complex_compliance = voigt_complex_comp + maxwell_complex_comp

    return complex_compliance


@tpy_vectorize(['complex128(float64, float64, float64, float64, float64)'])
def andrade(compliance, viscosity, frequency, alpha, zeta):
    """ Andrade Tides


    !TPY_args const: alpha, zeta
    """

    maxwell_complex_comp = maxwell(compliance, viscosity, frequency)

    andrade_term = compliance * viscosity * frequency * zeta
    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)
    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term
    imag_j[np.abs(andrade_term) <= float_eps] = 0.
    andrade_complex_comp = real_j + 1.0j * imag_j

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance


    if np.abs(frequency) <= float_eps:
        imag_j = 0.
    else:
        imag_j = -np.sin(alpha * np.pi / 2.) * const_term


    andrade_complex_comp = real_j + 1.0j * imag_j

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance


@tpy_vectorize(['complex128(float64, float64, float64, float64, float64, float64, float64)'])
def sundberg(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta):
    """ Sundberg-Cooper Tides

    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta

    """
    voigt_complex_comp = voigt(compliance, viscosity, frequency, voigt_compliance_offset, voigt_viscosity_offset)
    andrade_complex_comp = andrade(compliance, viscosity, frequency, alpha, zeta)

    complex_compliance = voigt_complex_comp + andrade_complex_comp

    return complex_compliance
