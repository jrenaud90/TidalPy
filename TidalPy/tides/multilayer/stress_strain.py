""" Calculation of stress, strain, and displacement based on tidal solutions and tidal potential

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

from typing import Tuple

import numpy as np

from .decompose import decompose
from ...utilities.performance import njit
from ...utilities.types import FloatArray

StressType = Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
StrainType = Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]


@njit(cacheable=True)
def calculate_displacements(
    tidal_potential: np.ndarray,
    tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
    tidal_solution_y: np.ndarray, colatitude: FloatArray
    ):
    """ Calculate tidal displacements using the tidal potential and its partial derivatives as well as the y-solution
    vector.

    Parameters
    ----------
    tidal_potential : FloatArray
        The tidal potential
    tidal_potential_partial_theta : FloatArray
        The partial derivative of the tidal potential with respect to the colatitude
    tidal_potential_partial_phi : FloatArray
        The partial derivative of the tidal potential with respect to the longitude
    tidal_solution_y : np.ndarray
        Vector of tidal solutions often denoted as a lower-case "y"
    colatitude : FloatArray
        Co-latitude where the calculations are conducted

    Returns
    -------
    radial_displacement : FloatArray
        The tidal displacement in the radial direction (r)
    polar_displacement : FloatArray
        The tidal displacement in the polar direction (theta)
    azimuthal_displacement : FloatArray
        The tidal displacement in the azimuthal direction (phi)
    """

    # TODO: add in summation over l and m?
    radial_displacement = tidal_solution_y[0, :] * tidal_potential
    polar_displacement = tidal_solution_y[2, :] * tidal_potential_partial_theta
    azimuthal_displacement = tidal_solution_y[2, :] * tidal_potential_partial_phi / np.sin(colatitude)

    return radial_displacement, polar_displacement, azimuthal_displacement


# @njit(cacheable=True)
def calculate_strain_stress_heating(
    tidal_potential: np.ndarray,
    tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
    tidal_potential_partial2_theta2: np.ndarray, tidal_potential_partial2_phi2: np.ndarray,
    tidal_potential_partial2_theta_phi: np.ndarray,
    tidal_solution_y: np.ndarray,
    colatitude: FloatArray,
    radius: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
    frequency: FloatArray, order_l: int = 2,
    ) -> Tuple[StrainType, StressType]:
    """ Calculate tidal strain tensor using the tidal potential and its partial derivatives as well as the y-solution
    vector.

    Parameters
    ----------
    tidal_potential : np.ndarray
        The tidal potential
    tidal_potential_partial_theta : np.ndarray
        The partial derivative of the tidal potential with respect to the colatitude
    tidal_potential_partial_phi : np.ndarray
        The partial derivative of the tidal potential with respect to the longitude
    tidal_potential_partial2_theta2 : np.ndarray
        The 2nd partial derivative of the tidal potential with respect to the colatitude
    tidal_potential_partial2_phi2 : np.ndarray
        The 2nd partial derivative of the tidal potential with respect to the longitude
    tidal_potential_partial2_theta_phi : np.ndarray
        The 2nd partial derivative of the tidal potential with respect to the colatitude and the longitude
    tidal_solution_y : np.ndarray
        Matrix [6 x N] of tidal solutions often denoted as a lower-case "y"
    colatitude : FloatArray
        Co-latitude where the calculations are conducted
    radius : np.ndarray
        Radius array [m]
    shear_moduli : np.ndarray
        Shear modulus as a function of radius [Pa]
    bulk_moduli : np.ndarray
        Bulk modulus as a function of radius [Pa]
    frequency : FloatArray
        Forcing frequency used to calculate tidal heating [rad s-1]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    strains : StrainType
        e_rr : np.ndarray
            Strain tensor component - radius radius
        e_thth : np.ndarray
            Strain tensor component - colatitude colatitude
        e_phph : np.ndarray
            Strain tensor component - longitude longitude
        e_rth : np.ndarray
            Strain tensor component - radius colatitude (appears twice in the symmetric tensor)
        e_rph : np.ndarray
            Strain tensor component - radius longitude (appears twice in the symmetric tensor)
        e_thph : np.ndarray
            Strain tensor component - colatitude longitude (appears twice in the symmetric tensor)
    stress : StressType
        s_rr : np.ndarray
            Stress tensor component - radius radius
        s_thth : np.ndarray
            Stress tensor component - colatitude colatitude
        s_phph : np.ndarray
            Stress tensor component - longitude longitude
        s_rth : np.ndarray
            Stress tensor component - radius colatitude (appears twice in the symmetric tensor)
        s_rph : np.ndarray
            Stress tensor component - radius longitude (appears twice in the symmetric tensor)
        s_thph : np.ndarray
            Stress tensor component - colatitude longitude (appears twice in the symmetric tensor)

    """

    # Optimizations
    sin_theta = np.sin(colatitude)
    cos_theta = np.cos(colatitude)
    cot_theta = cos_theta / sin_theta

    # Shortcuts
    r_inverse = (1. / radius)
    lame = bulk_moduli - (2. / 3.) * shear_moduli
    y1 = tidal_solution_y[0]
    y2 = tidal_solution_y[1]
    y3 = tidal_solution_y[2]
    y4 = tidal_solution_y[3]
    dy1_dr = (1. / (lame + 2. * shear_moduli)) * (y2 - (lame * r_inverse) * (2. * y1 - order_l * (order_l + 1.) * y3))

    # Strain Terms - TB05 Eq. 10; B13 Eq. 8 & 9
    # Build matrix for strain calculations
    # Strain and stress will have a shape of [6x strain/stress terms, radius_N, long_N, lat_N, time_N]
    strains = np.zeros((6, *y1.shape, *tidal_potential.shape), dtype=np.complex128)
    # \epsilon_{rr}
    strains[0, :, :, :, :] = dy1_dr * tidal_potential
    # \epsilon_{\theta\theta}
    strains[1, :, :, :, :] = r_inverse * (y3 * tidal_potential_partial2_theta2 + y1 * tidal_potential)
    # \epsilon_{\phi\phi}
    # There is a typo in Tobie+2005 with in both the \theta,\phi and \phi\phi components of the strain tensor,
    #    as pointed out in Kervazo et al (2021; A&A) Appendix D
    strains[2, :, :, :, :] = r_inverse * (y1 * tidal_potential +
                                          y3 * ((1. / sin_theta**2) * tidal_potential_partial2_phi2 +
                                                cot_theta * tidal_potential_partial_theta))
    # \epsilon_{r\theta}
    strains[3, :, :, :, :] = y4 * tidal_potential_partial_theta / shear_moduli
    # \epsilon_{r\phi}
    strains[4, :, :, :, :] = y4 * tidal_potential_partial_phi / (shear_moduli * sin_theta)
    # \epsilon_{\theta\phi}
    strains[5, :, :, :, :] = 2. * r_inverse * (y3 / sin_theta) * \
                            (tidal_potential_partial2_theta_phi - cot_theta * tidal_potential_partial_phi)

    # Calculate stress assuming isotropic media (Takeuchi & Saito (1972))
    # First calculate the trace of strain \epsilon_{rr} + \epsilon_{\theta\theta} + \epsilon_{\phi\phi}
    strain_trace = strains[0, :, :, :, :] + strains[1, :, :, :, :] + strains[2, :, :, :, :]

    # Constitutive equation: \sigma_{ij} = 2\bar{\mu} \epsilon_{ij} + \bar{\lambda}\epsilon_{kk}\delta_{ij}
    stresses = 2. * shear_moduli * strains
    # \epsilon_{rr}
    stresses[0, :, :, :, :] += lame * strain_trace
    # \epsilon_{\theta\theta}
    stresses[1, :, :, :, :] += lame * strain_trace
    # \epsilon_{\phi\phi}
    stresses[2, :, :, :, :] += lame * strain_trace

    # # Calculate Tidal Heating
    # volumetric_heating = (frequency / 2.) * (
    #         np.imag(s_rr) * np.real(e_rr) - np.real(s_rr) * np.imag(e_rr) +
    #         np.imag(s_thth) * np.real(e_thth) - np.real(s_thth) * np.imag(e_thth) +
    #         np.imag(s_phph) * np.real(e_phph) - np.real(s_phph) * np.imag(e_phph) +
    #         2. * (np.imag(s_rth) * np.real(e_rth) - np.real(s_rth) * np.imag(e_rth)) +
    #         2. * (np.imag(s_rph) * np.real(e_rph) - np.real(s_rph) * np.imag(e_rph)) +
    #         2. * (np.imag(s_thph) * np.real(e_thph) - np.real(s_thph) * np.imag(e_thph))
    # )

    # Uncomment to compare with Tobie method.
    # volumetric_heating, _ = decompose(
    #     tidal_solution_y, tidal_solution_y_derivative, radius, np.zeros_like(radius),
    #     shear_moduli, bulk_moduli, order_l=order_l
    #     )

    return strains, stresses
