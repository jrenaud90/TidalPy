""" Calculation of stress and strain based on tidal solutions and tidal potential.

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

from typing import Tuple

import numpy as np
from numba import prange

from ...utilities.performance import njit
from ...utilities.types import FloatArray

StressType = np.ndarray
StrainType = np.ndarray


@njit(cacheable=True, parallel=True)
def calculate_strain_stress(
    tidal_potential: np.ndarray,
    tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
    tidal_potential_partial2_theta2: np.ndarray, tidal_potential_partial2_phi2: np.ndarray,
    tidal_potential_partial2_theta_phi: np.ndarray,
    tidal_solution_y: np.ndarray,
    colatitude: FloatArray,
    radius: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
    frequency: FloatArray, order_l: int = 2
    ) -> Tuple[StrainType, StressType]:
    """ Calculate tidal strain tensor using the tidal potential and its partial derivatives as well as the y-solution
    vector.

    OPT: This function has actually been fairly optimized for performance since it is doing calculations on arrays that
        are n x m x k x p and can easily reach into the millions of elements. It also must be performed for each
        tidal mode which could be dozens per time step. It is pretty efficient right now but could always use a
        second pair of eyes due to its importance.

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
    y4_shear = y4 / shear_moduli
    dy1_dr = (1. / (lame + 2. * shear_moduli)) * (y2 - (lame * r_inverse) * (2. * y1 - order_l * (order_l + 1.) * y3))
    y3_r = y3 * r_inverse
    y1_r = y1 * r_inverse
    radius_n = y1.shape[0]
    stress_strain_shape = (6, radius_n, *tidal_potential.shape)

    # Build optimization terms
    s2_t1 = (1. / sin_theta**2) * tidal_potential_partial2_phi2 + cot_theta * tidal_potential_partial_theta
    s4_t0 = tidal_potential_partial_phi / sin_theta
    s5_t0 = 2. * (tidal_potential_partial2_theta_phi - cot_theta * tidal_potential_partial_phi) / sin_theta

    # Strain Terms - TB05 Eq. 10; B13 Eq. 8 & 9
    # Build matrix for strain calculations
    # Strain and stress will have a shape of [6x strain/stress terms, radius_N, long_N, lat_N, time_N]
    strains = np.empty(stress_strain_shape, dtype=np.complex128)
    stresses = np.empty(stress_strain_shape, dtype=np.complex128)

    for ri in prange(radius_n):
        # Pull out radius dependent parameters
        dy1dr_, y1_r_, y3_r_, y4_shear_, shear_, lame_ = \
            dy1_dr[ri], y1_r[ri], y3_r[ri], y4_shear[ri], shear_moduli[ri], lame[ri]

        # Any optimizations
        y1_r_tidal_potential = y1_r_ * tidal_potential

        # \epsilon_{rr}
        strains[0, ri, :, :, :] = dy1dr_ * tidal_potential
        # \epsilon_{\theta\theta}
        strains[1, ri, :, :, :] = y3_r_ * tidal_potential_partial2_theta2 + y1_r_tidal_potential
        # \epsilon_{\phi\phi}
        # There is a typo in Tobie+2005 with in both the \theta,\phi and \phi\phi components of the strain tensor,
        #    as pointed out in Kervazo et al (2021; A&A) Appendix D
        strains[2, ri, :, :, :] = y1_r_tidal_potential + y3_r_ * s2_t1

        # TODO: "Note that the non-diagonal strains (eh/, e/r, erh)in Takeuchi and Saito (1972)
        #  must be multiplied by 1/2 if they are to be the components of the strain tensor" Is this right?
        #  adding in the 1/2s now...
        # \epsilon_{r\theta}
        strains[3, ri, :, :, :] = y4_shear_ * tidal_potential_partial_theta / 2.
        # \epsilon_{r\phi}
        strains[4, ri, :, :, :] = y4_shear_ * s4_t0 / 2.
        # \epsilon_{\theta\phi}
        strains[5, ri, :, :, :] = y3_r_ * s5_t0 / 2.

        # Calculate stress assuming isotropic media (Takeuchi & Saito (1972))
        # First calculate the trace of strain \epsilon_{rr} + \epsilon_{\theta\theta} + \epsilon_{\phi\phi}
        strain_trace = strains[0, ri, :, :, :] + strains[1, ri, :, :, :] + strains[2, ri, :, :, :]
        strain_trace_lame = lame_ * strain_trace

        # Constitutive equation: \sigma_{ij} = 2\bar{\mu} \epsilon_{ij} + \bar{\lambda}\epsilon_{kk}\delta_{ij}
        stresses[:, ri, :, :, :] = (2. * shear_) * strains[:, ri, :, :, :]

        # \epsilon_{rr}
        stresses[0, ri, :, :, :] += strain_trace_lame
        # \epsilon_{\theta\theta}
        stresses[1, ri, :, :, :] += strain_trace_lame
        # \epsilon_{\phi\phi}
        stresses[2, ri, :, :, :] += strain_trace_lame

    # Calculate Tidal Heating
    # TODO: Does this frequency imply orbit average? There is still a time dependence...
    #   Do I take the orbit average here over this particular frequency cycle's period?
    # volumetric_heating = (frequency / 2.) * \
    #     np.sum(np.imag(stresses[:3]) * np.real(strains[:3]) - np.real(stresses[:3]) * np.imag(strains[:3]) +
    #            # 2 here is to double count cross terms since it is a symmetric matrix
    #            2. * (np.imag(stresses[3:]) * np.real(strains[3:]) - np.real(stresses[3:]) * np.imag(strains[3:])),
    #            axis=0)

    # Uncomment to compare with Tobie method.
    # volumetric_heating, _ = decompose(
    #     tidal_solution_y, tidal_solution_y_derivative, radius, np.zeros_like(radius),
    #     shear_moduli, bulk_moduli, order_l=order_l
    #     )

    return strains, stresses
