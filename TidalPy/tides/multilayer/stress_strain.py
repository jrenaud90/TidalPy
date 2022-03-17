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


@njit(cacheable=True)
def calculate_strain_stress_heating(
    tidal_potential: np.ndarray,
    tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
    tidal_potential_partial2_theta2: np.ndarray, tidal_potential_partial2_phi2: np.ndarray,
    tidal_potential_partial2_theta_phi: np.ndarray,
    tidal_solution_y: np.ndarray, tidal_solution_y_derivative: np.ndarray,
    colatitude: FloatArray,
    radius: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
    frequency: FloatArray, order_l: int = 2,
    ) -> Tuple[StrainType, StressType, np.ndarray]:
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
    tidal_solution_y_derivative: np.ndarray
        Matrix [6 x N] of the derivative of the tidal solutions with respect to radius (dy/dr).
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
    volumetric_heating : np.ndarray
        Volumetric heating rate [W]

    """

    # Optimizations
    sin_theta = np.sin(colatitude)
    cos_theta = np.cos(colatitude)
    cot_theta = cos_theta / sin_theta

    # Shortcuts
    dy1_dr = tidal_solution_y_derivative[0]
    y1 = tidal_solution_y[0]
    y2 = tidal_solution_y[1]
    y3 = tidal_solution_y[2]
    y4 = tidal_solution_y[3]

    # Build matrix for stress and strain calculations

    # Strain Terms - TB05 Eq. 10; B13 Eq. 8 & 9
    strain = np.zeros((6, *tidal_potential.shape), dtype=np.complex128)
    # \epsilon_{rr}
    strain[0, :, :, :, :] = dy1_dr * tidal_potential
    # \epsilon_{\theta\theta}
    strain[1, :, :, :, :] = (1. / radius) * (y3 * tidal_potential_partial2_theta2 + y1 * tidal_potential)
    # \epsilon_{\phi\phi}
    # There is a typo in Tobie+2005 with in both the \theta,\phi and \phi\phi components of the strain tensor,
    #    as pointed out in Kervazo et al (2021; A&A) Appendix D
    strain[2, :, :, :, :] = (1. / radius) * (y1 * tidal_potential +
                                             (y3 / sin_theta**2) * tidal_potential_partial2_phi2 +
                                             y3 * cot_theta * tidal_potential_partial_theta)
    # The (1/2) in the off-diagonal terms are due to these components appearing twice in the tensor
    # \epsilon_{r\theta}
    strain[3, :, :, :, :] = (1. / 2.) * y4 * tidal_potential_partial_theta / shear_moduli
    # \epsilon_{r\phi}
    strain[4, :, :, :, :] = (1. / 2.) * y4 * tidal_potential_partial_phi / (shear_moduli * sin_theta)
    # \epsilon_{\theta\phi}
    strain[5, :, :, :, :] = (1. / 2.) * (2. / radius) * (y3 / sin_theta) * \
                            (tidal_potential_partial2_theta_phi - cot_theta * tidal_potential_partial_phi)


    # Calculate stress using Kervazo et al (2021; A&A) Appendix D Eqs D.7-D.12
    stress = np.zeros((6, *tidal_potential.shape), dtype=np.complex128)
    strength_term_1 = (bulk_moduli - (2. / 3.) * shear_moduli)
    strength_term_2 = (bulk_moduli + (4. / 3.) * shear_moduli)
    theta_phi_term = strength_term_1 * dy1_dr + \
        (strength_term_2 / radius) * (2. * y1 - order_l * (order_l + 1.) * y3) - \
        2. * shear_moduli * y1 / radius

    # 2ue_ij + (K - 2/3u)*e_kk delta_ij

    2ue_rr + (K - 2/3u) * (e_rr + e_thth + e_phph)

    2u * dy1_dr * tidal_potential


    # \sigma_{rr}
    stress[0, :, :, :, :] = y2 * tidal_potential
    # \sigma_{\theta\theta}
    stress[1, :, :, :, :] = tidal_potential * (
            strength_term_1 * dy1_dr + (strength_term_2 / radius) * (2. * y1 - order_l * (order_l + 1.) * y3) -
            2. * shear_moduli * y1 / radius) - \
             (2. * shear_moduli * y3 / radius) * (cot_theta * tidal_potential_partial_theta +
                                                  (1. / sin_theta**2) * tidal_potential_partial2_phi2)
    # \epsilon_{\phi\phi}
    stress[2, :, :, :, :] = tidal_potential * (
            strength_term_1 * dy1_dr + (strength_term_2 / radius) * (2. * y1 - order_l * (order_l + 1.) * y3) -
                                2. * shear_moduli * y1 / radius) - \
             (2. * shear_moduli * y3 / radius) * tidal_potential_partial2_theta2
    s_thph = (2. * shear_moduli * y3 / (radius * sin_theta)) * (tidal_potential_partial2_theta_phi -
                                                                cot_theta * tidal_potential_partial_phi)
    s_rth = y4 * tidal_potential_partial_theta
    s_rph = (y4 / sin_theta) * tidal_potential_partial_phi

    # The (1/2) in the off-diagonal terms are due to these components appearing twice in the tensor
    s_rph *= (1. / 2.)
    s_rth *= (1. / 2.)
    s_thph *= (1. / 2.)

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

    return strains, stresses, volumetric_heating
