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
import warnings

import numpy as np

from ...utilities.performance import find_factorial, njit
from ...utilities.types import ComplexArray, FloatArray, float_eps, float_lognat_max

warnings.warn('Deprecation Warning: the non-cythonized TidalPy.rheology.complex_compliance.compliance_models will be removed in a future release of TidalPy. Please use TidalPy.rheology.models instead. Please report any differences noted so that they can be addressed before the future release of TidalPy.', DeprecationWarning)


# OPT: # TODO: @vectorize(['complex128(float64, float64, float64)'],nopython=True) seems to do everything we need for
#    these functions and then we can return to the if/else version for frequency. speeds are better when dealing wit
#    large arrays. The downside, and this is a big downside, is the speed is much slower for all float (non-arrays)
#    which are important for differential equations.

@njit(cacheable=False)
def off(frequency: 'FloatArray', compliance: 'FloatArray', viscosity: FloatArray) -> 'ComplexArray':
    """ Calculates the complex compliance utilizing the model: Off

    !TPY_args live: self.compliance, self.viscosity

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    shape = viscosity + compliance + frequency

    real_j = compliance

    complex_compliance = real_j + \
                         (shape <= float_eps) * 0.0j

    return complex_compliance


@njit(cacheable=False)
def fixed_q(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    planet_beta: float = 3.443e11, quality_factor: float = 10.
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Fixed-Q

    !TPY_args live: self.compliance, self.viscosity, self.beta, self.quality_factor

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]
    planet_beta : FloatArray
        Planet or Layer: Radius * Density * Gravity Accl. [kg m-1 s-2]
    quality_factor : FloatArray
        Planet or Layer's Quality factor (k_2 / Q)

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    shape = 0. * (viscosity + compliance + frequency)

    real_j = -19. / (2. * planet_beta)
    imag_j = -quality_factor * (19. / 4.) * (2. * planet_beta * compliance + 19.) / (compliance * planet_beta**2) \
             * (1. + 0. * frequency)

    complex_compliance = real_j + \
                         ((np.abs(frequency) + shape) <= float_eps) * 0.0j + \
                         ((np.abs(frequency) + shape) > float_eps) * 1.0j * imag_j

    return complex_compliance

@njit(cacheable=False)
def newton(frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray') -> 'ComplexArray':
    """ Calculates the complex compliance utilizing the model: Newton

    !TPY_args live: self.compliance, self.viscosity

    Notes
    -----

    References
    ----------

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    shape = 0. * (viscosity + compliance + frequency)

    denominator = (viscosity * frequency)
    denominator = (np.abs(denominator) <= float_eps) * 1.0e-100 + \
                  (np.abs(denominator) > float_eps) * denominator

    real_j = 0.
    imag_j = -1.0 / denominator

    complex_compliance = real_j + \
                         ((np.abs(frequency) + shape) <= float_eps) * 0.0j + \
                         ((np.abs(frequency) + shape) > float_eps) * 1.0j * imag_j

    return complex_compliance

@njit(cacheable=False)
def elastic(frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray') -> 'ComplexArray':
    """ Calculates the complex compliance utilizing the model: Elastic

    !TPY_args live: self.compliance, self.viscosity

    Notes
    -----

    References
    ----------

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    shape = 0. * (viscosity + compliance + frequency)

    denominator = (viscosity * frequency)
    denominator = (np.abs(denominator) <= float_eps) * 1.0e-100 + \
                  (np.abs(denominator) > float_eps) * denominator

    real_j = compliance
    imag_j = 0.

    complex_compliance = real_j + \
                         ((np.abs(frequency) + shape) <= float_eps) * 0.0j + \
                         ((np.abs(frequency) + shape) > float_eps) * 1.0j * imag_j

    return complex_compliance


@njit(cacheable=False)
def maxwell(frequency: 'FloatArray', compliance: 'FloatArray', viscosity: FloatArray) -> 'ComplexArray':
    """ Calculates the complex compliance utilizing the model: Maxwell

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
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    shape = 0. * (viscosity + compliance + frequency)

    denominator = (viscosity * frequency)
    denominator = (np.abs(denominator) <= float_eps) * 1.0e-100 + \
                  (np.abs(denominator) > float_eps) * denominator

    real_j = compliance
    imag_j = -1.0 / denominator

    complex_compliance = real_j + \
                         ((np.abs(frequency) + shape) <= float_eps) * 0.0j + \
                         ((np.abs(frequency) + shape) > float_eps) * 1.0j * imag_j

    return complex_compliance


@njit(cacheable=False)
def voigt(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Voigt-Kelvin

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
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    shape = 0. * (viscosity + compliance + frequency)

    voigt_comp = voigt_compliance_offset * compliance
    voigt_visc = voigt_viscosity_offset * viscosity

    denominator = (voigt_comp * voigt_visc * frequency)**2 + 1.
    denominator = (np.abs(denominator) <= float_eps) * 1.0e-100 + \
                  (np.abs(denominator) > float_eps) * denominator

    real_j = voigt_comp / denominator
    imag_j = -voigt_comp**2 * voigt_visc * frequency / denominator

    complex_compliance = real_j + \
                         ((np.abs(frequency) + shape) <= float_eps) * 0.0j + \
                         ((np.abs(frequency) + shape) > float_eps) * 1.0j * imag_j

    return complex_compliance


@njit(cacheable=False)
def burgers(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Burgers

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
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    voigt_complex_comp = voigt(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)
    maxwell_complex_comp = maxwell(frequency, compliance, viscosity)

    complex_compliance = voigt_complex_comp + maxwell_complex_comp

    return complex_compliance


@njit(cacheable=False)
def andrade(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    alpha: float = 0.3, zeta: float = 1.
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Andrade

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: alpha, zeta

    Notes
    -----
    The Andrade rheology is partially constructed from the Maxwell rheology. This is further modified by a term that,
    in the time-domain, is proportional to t^{\alpha}. In the Fourier domain this translates to a frequency dependence
    ~\omega^{-\alpha}. This model was originally developed for the stress-strain relationship in metals, but has been
    found to model planetary materials as well.

    This version of the model will not transition into a Maxwell-like rheology at very low frequencies.

    References
    ----------
    - Gribb and Cooper (1998), JGR, DOI: 10.1029/98JB02786
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    andrade_term = compliance * viscosity * frequency * zeta
    andrade_term = (np.abs(andrade_term) <= float_eps) * 1.0e-100 + \
                   (np.abs(andrade_term) > float_eps) * andrade_term

    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)

    shape = 0. * (andrade_term + alpha)

    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term

    # If the frequency is at zero then the real value goes to +infinity. Set to large value instead. Imag -> 0.
    andrade_complex_comp = ((np.abs(frequency) + shape) <= float_eps) * (1.0e100 + 0.0j) + \
                           ((np.abs(frequency) + shape) > float_eps) * (real_j + 1.0j * imag_j)

    maxwell_complex_comp = maxwell(frequency, compliance, viscosity)

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance


@njit(cacheable=False)
def andrade_freq(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    alpha: float = 0.3, zeta: float = 1., critical_freq: float = 7.27221e-7, critical_freq_falloff: float = 30
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Andrade with a frequency-dependent zeta

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: alpha, zeta, critical_freq, critical_freq_falloff

    Notes
    -----
    The Andrade rheology is partially constructed from the Maxwell rheology. This is further modified by a term that,
    in the time-domain, is proportional to t^{\alpha}. In the Fourier domain this translates to a frequency dependence
    ~\omega^{-\alpha}. This model was originally developed for the stress-strain relationship in metals, but has been
    found to model planetary materials as well.

    This version of the model will transition into a Maxwell-like rheology at very low frequencies.

    References
    ----------
    - Karato and Spetzler (1990), RGeo, DOI: 10.1029/RG028i004p00399
    - Gribb and Cooper (1998), JGR, DOI: 10.1029/98JB02786
    - Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter equal to the Andrade characteristic timescale / Maxwell
    critical_freq : float = 7.27221e-7
        For forcing frequencies smaller than the critical frequency, the Andrade component converges to
        the Maxwell component [rads s-1].
        Default is 2 * pi / 100 days.
    critical_freq_falloff : float = 30
        Determines how quickly the Andrade rheology turns to Maxwell falls off after dropping below the critical freq.

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    # Update Zeta based on an additional frequency dependence.
    # We will say that zeta can not get above 1e100 (it would be Maxwell like well before that)
    freq_ratio = np.abs(frequency / critical_freq)
    exponent = -critical_freq_falloff * (freq_ratio - 1.)
    exponent = (exponent >= 100.) * 100. + \
               (exponent <= 0.) * 0. + \
               (exponent > 0.) * (exponent < 100.) * exponent

    updated_zeta = zeta * np.exp(exponent)

    # Continue on with regular Andrade calculation
    andrade_term = compliance * viscosity * frequency * updated_zeta
    andrade_term = (np.abs(andrade_term) <= float_eps) * 1.0e-100 + \
                   (np.abs(andrade_term) > float_eps) * andrade_term

    const_term = compliance * andrade_term**(-alpha) * find_factorial(alpha)

    shape = 0. * (andrade_term + alpha)

    real_j = np.cos(alpha * np.pi / 2.) * const_term
    imag_j = -np.sin(alpha * np.pi / 2.) * const_term

    # If the frequency is at zero then the real value goes to +infinity. Set to large value instead. Imag -> 0.
    andrade_complex_comp = ((np.abs(frequency) + shape) <= float_eps) * (1.0e100 + 0.0j) + \
                           ((np.abs(frequency) + shape) > float_eps) * (real_j + 1.0j * imag_j)

    maxwell_complex_comp = maxwell(frequency, compliance, viscosity)

    complex_compliance = maxwell_complex_comp + andrade_complex_comp

    return complex_compliance


@njit(cacheable=False)
def sundberg(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02,
    alpha: float = 0.3, zeta: float = 1.
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Sundberg-Cooper

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta

    Notes
    -----
    The Sundberg-Cooper rheology is a linear sum of the Andrade and Burgers rheologies. However, even though its
    Parameters share the same symbol and names, they may differ from those used for either of its composite model.

    This version of the model will not transition into a Burgers-like rheology at very low frequencies.

    References
    ----------
    - Sundberg and Cooper (2010), Philo. Mag., DOI: 10.1080/14786431003746656
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
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
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    andrade_complex_comp = \
        andrade(frequency, compliance, viscosity, alpha, zeta)

    voigt_complex_comp = \
        voigt(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)

    complex_compliance = voigt_complex_comp + andrade_complex_comp

    return complex_compliance


@njit(cacheable=False)
def sundberg_freq(
    frequency: 'FloatArray', compliance: 'FloatArray', viscosity: 'FloatArray',
    voigt_compliance_offset: float = 0.2, voigt_viscosity_offset: float = 0.02,
    alpha: float = 0.3, zeta: float = 1.,  critical_freq: float = 7.27221e-7, critical_freq_falloff: float = 30
    ) -> ComplexArray:
    """ Calculates the complex compliance utilizing the model: Sundberg-Cooper with a frequency-dependent zeta

    !TPY_args live: self.compliance, self.viscosity
    !TPY_args const: voigt_compliance_offset, voigt_viscosity_offset, alpha, zeta, critical_freq, critical_freq_falloff

    Notes
    -----
    The Sundberg-Cooper rheology is a linear sum of the Andrade and Burgers rheologies. However, even though its
    Parameters share the same symbol and names, they may differ from those used for either of its composite model.

    This version of the model will transition into a Burgers-like rheology at very low frequencies.

    References
    ----------
    - Karato and Spetzler (1990), RGeo, DOI: 10.1029/RG028i004p00399
    - Sundberg and Cooper (2010), Philo. Mag., DOI: 10.1080/14786431003746656
    - Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784

    Parameters
    ----------
    frequency : FloatArray
        Tidal forcing frequency [rads s-1]
        Note that a world may experience multiple tidal frequencies for NSR tides, if the eccentricity or obliquity
        is large, or for a tidal harmonic integer l > 2.
    compliance : FloatArray
        Layer or Planet's compliance (inverse of shear modulus) [Pa-1]
    viscosity : FloatArray
        Layer or Planet's effective viscosity [Pa s]
    voigt_compliance_offset : float
        Voigt component's compliance offset eta_voigt = voigt_compliance_offset * compliance
    voigt_viscosity_offset : float
        Voigt component's viscosity offset eta_voigt = voigt_viscosity_offset * viscosity
    alpha : float
        Andrade exponent parameter
    zeta : float
        Andrade timescale parameter equal to the Andrade characteristic timescale / Maxwell
    critical_freq : float = 7.27221e-7
        For forcing frequencies smaller than the critical frequency, the Andrade component converges to
        the Maxwell component [rads s-1].
        Default is 2 * pi / 100 days.
    critical_freq_falloff : float = 30
        Determines how quickly the Andrade rheology turns to Maxwell falls off after dropping below the critical freq.

    Returns
    -------
    complex_compliance : ComplexArray
        Complex compliance (complex number) [Pa-1]
    """

    andrade_freq_complex_comp = \
        andrade_freq(frequency, compliance, viscosity, alpha, zeta, critical_freq, critical_freq_falloff)

    voigt_complex_comp = \
        voigt(frequency, compliance, viscosity, voigt_compliance_offset, voigt_viscosity_offset)

    complex_compliance = voigt_complex_comp + andrade_freq_complex_comp

    return complex_compliance

# Put New Models Below Here!
