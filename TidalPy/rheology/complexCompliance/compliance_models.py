""" Tides Compliance Functions

How To Implement a new compliance model:
    * Make a new function here with a function doc-string that follows the format of the other pre-built
        complex compliance functions found here.

        - All complex compliance functions *require* at least three arguments: compliance (non-complex), viscosity,
            and frequency. These should be the first arguments and left in this order.

        - Additional arguments can be supplied by adding them after these required arguments but they must be
            included in the function's doc-string using either (or both) the "!TPY_args live:" or "!TPY_args const:"
            prefix. The arguments should be listed in the doc-string in the same order they are in the argument field.

    * Any new arguments that are simple constants can then be added to the rheology section of a planet's layer
        configuration. These should have a "!TPY_args const:" prefix. Non-constants require "!TPY_args live:" prefix.

    * You will also need to make a function that works with arrays. If the above function works fine with both
        floats and arrays then you don't need to do anything. Otherwise you should make a float safe version
        called "name" and an array safe one called "name_array"
"""

import numpy as np

from ...utilities.performance.numba import find_factorial, njit
from ...utilities.types import float_eps, float_lognat_max


@njit
def off(frequency: float, compliance: float, viscosity: float) -> complex:
    """ Calculates the complex compliance utilizing the model: Off - NonArray Only

    !TPY_args live: self.compliance, self.viscosity

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    real_j = compliance
    imag_j = 0.

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def off_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Off - Arrays only

    !TPY_args live: self.compliance, self.viscosity

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    real_j = compliance
    imag_j = np.zeros_like(frequency + viscosity + compliance)

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def fixed_q(frequency: float, compliance: float, viscosity: float,
            planet_beta: float = 3.443e11, quality_factor: float = 10.) -> complex:
    """ Calculates the complex compliance utilizing the model: Fixed-Q - NonArray Only

    !TPY_args live: self.compliance, self.viscosity, self.beta, self.quality_factor

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    planet_beta : float
        Planet or Layer: Radius * Density * Gravity Accl. [kg m-1 s-2]
    quality_factor : float
        Planet or Layer's Quality factor (k_2 / Q)

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    real_j = -19. / (2. * planet_beta)

    if abs(frequency) <= float_eps:
        imag_j = 0.
    else:
        imag_j = -quality_factor * (19. / 4.) * (2. * planet_beta * compliance + 19.) / (compliance * planet_beta**2) \
                 * (1. + 0. * frequency)

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def fixed_q_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
                  planet_beta: float = 3.443e11, quality_factor: float = 10.) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Fixed-Q - Arrays only

    !TPY_args live: self.compliance, self.viscosity, self.beta, self.quality_factor

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    planet_beta : float
        Planet or Layer: Radius * Density * Gravity Accl. [kg m-1 s-2]
    quality_factor : float
        Planet or Layer's Quality factor (k_2 / Q)

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    real_j = -19. / (2. * planet_beta)
    # The (1 + 0*w) hack is there to ensure that imag_j shares shape with the larger of the two: compliance or frequency
    imag_j = -quality_factor * (19. / 4.) * (2. * planet_beta * compliance + 19.) / (compliance * planet_beta**2) \
        * (1. + 0. * frequency)
    imag_j[np.abs(frequency * np.ones_like(compliance * viscosity)) <= float_eps] = 0.

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance


@njit
def maxwell(frequency: float, compliance: float, viscosity: float) -> complex:
    """ Calculates the complex compliance utilizing the model: Maxwell - NonArray Only

    !TPY_args live: self.compliance, self.viscosity

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    real_j = compliance
    denominator = (viscosity * frequency)

    if abs(frequency) <= float_eps:
        imag_j = 0.
    else:
        imag_j = -1.0 / denominator

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def maxwell_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Maxwell - Array only

    !TPY_args live: self.compliance, self.viscosity

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    real_j = compliance
    denominator = (viscosity * frequency)
    imag_j = -1.0 / denominator
    imag_j[np.abs(denominator) <= float_eps] = 0.

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def voigt(frequency: float, compliance: float, viscosity: float,
          voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02) -> complex:
    """ Calculates the complex compliance utilizing the model: Voigt-Kelvin - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    voigt_comp = voigt_compliance_offset * compliance
    voigt_visc = voigt_viscosity_offset * viscosity

    denominator = (voigt_comp * voigt_visc * frequency)**2 + 1.
    real_j = voigt_comp / denominator
    if abs(frequency) <= float_eps:
        imag_j = 0.
    else:
        imag_j = -voigt_comp**2 * voigt_visc * frequency / denominator

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def voigt_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
                voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Voigt-Kelvin - Arrays Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    voigt_comp = voigt_compliance_offset * compliance
    voigt_visc = voigt_viscosity_offset * viscosity

    denominator = (voigt_comp * voigt_visc * frequency)**2 + 1.
    real_j = voigt_comp / denominator
    imag_j = -voigt_comp**2 * voigt_visc * frequency / denominator
    imag_j[np.abs(frequency) <= float_eps] = 0.

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def burgers(frequency: float, compliance: float, viscosity: float,
            voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02) -> complex:
    """ Calculates the complex compliance utilizing the model: Burgers - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    voigt_complex_comp = voigt(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)
    maxwell_complex_comp = maxwell(frequency, compliance, viscosity)

    complex_compliance = voigt_complex_comp + maxwell_complex_comp

    return complex_compliance

@njit
def burgers_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
                  voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Burgers - Arrays Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    voigt_complex_comp = voigt_array(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)
    maxwell_complex_comp = maxwell_array(frequency, compliance, viscosity)

    complex_compliance = voigt_complex_comp + maxwell_complex_comp

    return complex_compliance

@njit
def andrade(frequency: float, compliance: float, viscosity: float,
            alpha: float = 0.3, zeta: float = 1.) -> complex:
    """ Calculates the complex compliance utilizing the model: Andrade - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: alpha, zeta

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    andrade_term = compliance * viscosity * frequency * zeta
    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)

    if abs(frequency) <= float_eps:
        real_j = 1.e100
        imag_j = 0.
    else:
        real_j =  np.cos(alpha * np.pi / 2.) * const_term
        imag_j = -np.sin(alpha * np.pi / 2.) * const_term

    andrade_complex_comp = real_j + 1.0j * imag_j
    maxwell_complex_comp = maxwell(frequency, compliance, viscosity)

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def andrade_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
                  alpha: float = 0.3, zeta: float = 1.) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Andrade - Arrays Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: alpha, zeta

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    andrade_term = compliance * viscosity * frequency * zeta
    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)
    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term

    real_j[np.abs(andrade_term) <= float_eps] = 1.e100
    imag_j[np.abs(andrade_term) <= float_eps] = 0.
    andrade_complex_comp = real_j + 1.0j * imag_j

    maxwell_complex_comp = maxwell_array(frequency, compliance, viscosity)

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def andrade_freq(frequency: float, compliance: float, viscosity: float,
                 alpha: float = 0.3, zeta: float = 1., critical_freq: float = 2.e-5) -> complex:
    """ Calculates the complex compliance utilizing the model: Andrade with a frequency-dependent zeta - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: alpha, zeta, critical_freq

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter
    critical_freq : float
        For forcing frequencies larger than the critical frequency, the Andrade component converges to
            the Maxwell component [rads s-1]

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    # Update Zeta based on an additional frequency dependence.
    freq_ratio = abs(critical_freq / frequency)
    if freq_ratio > float_lognat_max:
        freq_ratio = float_lognat_max

    zeta *= np.exp(freq_ratio)

    andrade_term = compliance * viscosity * frequency * zeta
    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)
    if abs(frequency) <= float_eps:
        real_j = 1.e100
        imag_j = 0.
    else:
        real_j = np.cos(alpha*np.pi/2.)*const_term
        imag_j = -np.sin(alpha*np.pi/2.)*const_term
    andrade_complex_comp = real_j + 1.0j * imag_j

    maxwell_complex_comp = maxwell(frequency, compliance, viscosity)

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def andrade_freq_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
                 alpha: float = 0.3, zeta: float = 1., critical_freq: float = 2.e-5) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Andrade with a frequency-dependent zeta - Arrays Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: alpha, zeta, critical_freq

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter
    critical_freq : float
        For forcing frequencies larger than the critical frequency, the Andrade component converges to
            the Maxwell component [rads s-1]

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    # Update Zeta based on an additional frequency dependence.
    freq_ratio = np.abs(critical_freq / frequency)
    freq_ratio[freq_ratio > float_lognat_max] = float_lognat_max
    zeta = zeta * np.exp(freq_ratio)

    andrade_term = compliance * viscosity * frequency * zeta
    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)
    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term

    real_j[np.abs(andrade_term) <= float_eps] = 1.e100
    imag_j[np.abs(andrade_term) <= float_eps] = 0.
    andrade_complex_comp = real_j + 1.0j * imag_j

    maxwell_complex_comp = maxwell_array(frequency, compliance, viscosity)

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def sundberg(frequency: float, compliance: float, viscosity: float,
             voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02,
             alpha: float = 0.3, zeta: float = 1.) -> complex:
    """ Calculates the complex compliance utilizing the model: Sundberg-Cooper - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    andrade_complex_comp = andrade(frequency, compliance, viscosity, alpha, zeta)
    voigt_complex_comp = voigt(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)

    complex_compliance = voigt_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def sundberg_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
             voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02,
             alpha: float = 0.3, zeta: float = 1.) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Sundberg-Cooper - Arrays Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    andrade_complex_comp = andrade_array(frequency, compliance, viscosity, alpha, zeta)
    voigt_complex_comp = voigt_array(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)

    complex_compliance = voigt_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def sundberg_freq(frequency: float, compliance: float, viscosity: float,
                  voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02,
                  alpha: float = 0.3, zeta: float = 1., critical_freq: float = 2.e-5) -> complex:
    """ Calculates the complex compliance utilizing the model: Sundberg-Cooper with a frequency-dependent zeta - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta, critical_freq

    Parameters
    ----------
    frequency : float
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : float
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : float
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter
    critical_freq : float
        For forcing frequencies larger than the critical frequency, the Andrade component converges to
            the Maxwell component [rads s-1]

    Returns
    -------
    complex_compliance : complex
        Complex compliance (complex number) [Pa-1]
    """

    andrade_complex_comp = andrade_freq(frequency, compliance, viscosity, alpha, zeta, critical_freq)
    voigt_complex_comp = voigt(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)

    complex_compliance = voigt_complex_comp + andrade_complex_comp

    return complex_compliance

@njit
def sundberg_freq_array(frequency: np.ndarray, compliance: np.ndarray, viscosity: np.ndarray,
                        voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02,
                        alpha: float = 0.3, zeta: float = 1., critical_freq: float = 2.e-5) -> np.ndarray:
    """ Calculates the complex compliance utilizing the model: Sundberg-Cooper with a frequency-dependent zeta - Arrays Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta, critical_freq

    Parameters
    ----------
    frequency : np.ndarray
        Planet's tidal frequency [rads s-1]
        Note that a planet may experience multiple tidal frequencies for NSR tides or Fourier Degree l>2
    compliance : np.ndarray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : np.ndarray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter
    critical_freq : float
        For forcing frequencies larger than the critical frequency, the Andrade component converges to
            the Maxwell component [rads s-1]

    Returns
    -------
    complex_compliance : np.ndarray
        Complex compliance (complex number) [Pa-1]
    """

    andrade_complex_comp = andrade_freq_array(frequency, compliance, viscosity, alpha, zeta, critical_freq)
    voigt_complex_comp = voigt_array(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)

    complex_compliance = voigt_complex_comp + andrade_complex_comp

    return complex_compliance

# Put New Models Below Here!