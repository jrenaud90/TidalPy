""" Complex Compliance Functions

This module provides functions to calculate the complex compliance for various rheological models. The complex
compliance is the inverse of the complex shear modulus and is used to estimate tidal dissipation. It depends on the
static (non-complex) shear modulus, viscosity, forcing frequency, and other material properties --- some of which may
depend upon temperature, melt-fraction, etc. The functional form of this dependence is determined by the material's
rheology.

Notes
-----
    How To Implement a new complex compliance (rheology) model:
        - Make a new function here with a function doc-string that follows the format of the other pre-built
          complex compliance functions found here.
            - All complex compliance functions *require* at least three arguments: frequency, compliance (non-complex),
              viscosity. These must be the first arguments and kept in this order.
            - Additional arguments can be supplied by adding them after these required arguments but they must be
              included in the function's doc-string using either (or both) the "!TPY_args live:" or "!TPY_args const:"
              prefix. The arguments should be listed in the doc-string in the same order they are in the argument field.
        - Any new arguments that are simple constants can then be added to the rheology section of a planet's layer
          configuration. These should have a "!TPY_args const:" prefix. Non-constants require "!TPY_args live:" prefix.
        - You will also need to make a function that works with arrays. If the above function works fine with both
          floats and arrays then you don't need to do anything. Otherwise you should make a float safe version
          called "name" and an array safe one called "name_array"

References
----------
Rheologies and the relationship between stress and strain have had literally hundreds of years of study. However, the
following references provide a good starting point for the currently implemented rheologies in TidalPy in the context
of tides.

    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000
        - Background information on the Maxwell, Voigt-Kelvin, and Burgers rheologies.
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
        - Details on how complex compliances are related to global Love numbers.
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784
        - Background information on the Andrade and Sundberg-Cooper rheologies.


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
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Maxwell rheology has been the traditional model used to estimate tidal dissipation in planets and moons.
    It has a characteristic peak in dissipation vs. forcing frequency. However, it has been shown that it
    underestimates dissipation at high frequencies.

    References
    ----------
    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Maxwell rheology has been the traditional model used to estimate tidal dissipation in planets and moons.
    It has a characteristic peak in dissipation vs. forcing frequency. However, it has been shown that it
    underestimates dissipation at high frequencies.

    References
    ----------
    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Voigt-Kelvin rheology is a non-realistic (except for very specific circumstances) that is generally more useful
    when comparing or building other models. It is characterized by an 'island' of dissipation in the shear modulus vs.
    viscosity phase space.

    References
    ----------
    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Voigt-Kelvin rheology is a non-realistic (except for very specific circumstances) that is generally more useful
    when comparing or building other models. It is characterized by an 'island' of dissipation in the shear modulus vs.
    viscosity phase space.

    References
    ----------
    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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
    imag_j[np.abs(denominator) <= float_eps] = 0.

    complex_compliance = real_j + 1.0j * imag_j

    return complex_compliance

@njit
def burgers(frequency: float, compliance: float, viscosity: float,
            voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02) -> complex:
    """ Calculates the complex compliance utilizing the model: Burgers - NonArray Only

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset

    Notes
    -----
    The Burgers rheology exhibits a secondary peak in the dissipation vs. frequency domain. This peak describes a
    secondary dissipation mechanism that is dominant at this forcing frequency (grain boundary sliding, dislocation
    diffusion, etc.). It can be constructed by a linear summation of the Voigt-Kelvin and Maxwell rheologies.

    References
    ----------
    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Burgers rheology exhibits a secondary peak in the dissipation vs. frequency domain. This peak describes a
    secondary dissipation mechanism that is dominant at this forcing frequency (grain boundary sliding, dislocation
    diffusion, etc.).

    References
    ----------
    - Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Andrade rheology is partially constructed from the Maxwell rheology. This is further modified by a term that,
    in the time-domain, is proportional to t^{\alpha}. In the Fourier domain this translates to a frequency dependence
    ~\omega^{-\alpha}. This model was originally developed for the stress-strain relationship in metals, but has been
    found to model planetary materials as well.

    References
    ----------
    - Gribb and Cooper (1998), JGR, DOI: 10.1029/98JB02786
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Andrade rheology is partially constructed from the Maxwell rheology. This is further modified by a term that,
    in the time-domain, is proportional to t^{\alpha}. In the Fourier domain this translates to a frequency dependence
    ~\omega^{-\alpha}. This model was originally developed for the stress-strain relationship in metals, but has been
    found to model planetary materials as well.

    References
    ----------
    - Gribb and Cooper (1998), JGR, DOI: 10.1029/98JB02786
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Andrade rheology is partially constructed from the Maxwell rheology. This is further modified by a term that,
    in the time-domain, is proportional to t^{\alpha}. In the Fourier domain this translates to a frequency dependence
    ~\omega^{-\alpha}. This model was originally developed for the stress-strain relationship in metals, but has been
    found to model planetary materials as well.

    References
    ----------
    - Gribb and Cooper (1998), JGR, DOI: 10.1029/98JB02786
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Andrade rheology is partially constructed from the Maxwell rheology. This is further modified by a term that,
    in the time-domain, is proportional to t^{\alpha}. In the Fourier domain this translates to a frequency dependence
    ~\omega^{-\alpha}. This model was originally developed for the stress-strain relationship in metals, but has been
    found to model planetary materials as well.

    References
    ----------
    - Gribb and Cooper (1998), JGR, DOI: 10.1029/98JB02786
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Sundberg-Cooper rheology is a linear sum of the Andrade and Burgers rheologies. However, even though its
    Parameters share the same symbol and names, they may differ from those used for either of its composite model.

    References
    ----------
    - Sundberg and Cooper (2010), Philo. Mag., DOI: 10.1080/14786431003746656
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Sundberg-Cooper rheology is a linear sum of the Andrade and Burgers rheologies. However, even though its
    Parameters share the same symbol and names, they may differ from those used for either of its composite model.

    References
    ----------
    - Sundberg and Cooper (2010), Philo. Mag., DOI: 10.1080/14786431003746656
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Sundberg-Cooper rheology is a linear sum of the Andrade and Burgers rheologies. However, even though its
    Parameters share the same symbol and names, they may differ from those used for either of its composite model.

    References
    ----------
    - Sundberg and Cooper (2010), Philo. Mag., DOI: 10.1080/14786431003746656
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : float
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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

    Notes
    -----
    The Sundberg-Cooper rheology is a linear sum of the Andrade and Burgers rheologies. However, even though its
    Parameters share the same symbol and names, they may differ from those used for either of its composite model.

    References
    ----------
    - Sundberg and Cooper (2010), Philo. Mag., DOI: 10.1080/14786431003746656
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : np.ndarray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
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