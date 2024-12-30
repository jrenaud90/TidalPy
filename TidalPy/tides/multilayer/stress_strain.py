""" Calculation of stress and strain based on tidal solutions and tidal potential.

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

from typing import Tuple, TYPE_CHECKING

import numpy as np

from TidalPy.utilities.performance import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray

StressType = np.ndarray
StrainType = np.ndarray


@njit(cacheable=True)
def calculate_strain_stress(
    tidal_potential: np.ndarray,
    tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
    tidal_potential_partial2_theta2: np.ndarray, tidal_potential_partial2_phi2: np.ndarray,
    tidal_potential_partial2_theta_phi: np.ndarray,
    tidal_solution_y: np.ndarray,
    longitude_array: np.ndarray, colatitude_array: np.ndarray, time_array: np.ndarray,
    radius_array: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
    frequency: 'FloatArray', degree_l: int = 2
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
    longitude_array : np.ndarray
        Array of longitudes where tides were calculated [radians]
    colatitude_array : np.ndarray
        Array of colatitudes where tides were calculated [radians]
    time_array : np.ndarray
        Array of points in time when tides were calculated [s]
    radius_array : np.ndarray
        Radius array [m]
    shear_moduli : np.ndarray
        Shear modulus as a function of radius [Pa]
    bulk_moduli : np.ndarray
        Bulk modulus as a function of radius [Pa]
    frequency : FloatArray
        Forcing frequency used to calculate tidal heating [rad s-1]
    degree_l : int = 2
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
    # Get size of arrays
    n_radius = len(radius_array)
    n_longitude = len(longitude_array)
    n_colatitude = len(colatitude_array)
    n_time = len(time_array)

    # Strain Terms - TB05 Eq. 10; B13 Eq. 8 & 9
    # Build matrix for strain calculations
    # Strain and stress will have a shape of [6x strain/stress terms, radius_N, long_N, lat_N, time_N]
    strains = np.empty((6, n_radius, n_longitude, n_colatitude, n_time), dtype=np.complex128)
    stresses = np.empty((6, n_radius, n_longitude, n_colatitude, n_time), dtype=np.complex128)

    for ri in range(n_radius):
        # Pull out radius dependent parameters
        radius = radius_array[ri]
        shear  = shear_moduli[ri]
        bulk   = bulk_moduli[ri]
        y1     = tidal_solution_y[0, ri]
        y2     = tidal_solution_y[1, ri]
        y3     = tidal_solution_y[2, ri]
        y4     = tidal_solution_y[3, ri]
        
        # Optimizations
        lame = bulk - (2. / 3.) * shear
        lame_2shear = (lame + 2. * shear)
        if shear == 0.0:
            y4_shear   = 0.0
        else:
            y4_shear   = y4 / shear

        if radius == 0.0:
            radius_inv = np.nan
            dy1_dr     = np.nan
            y3_r       = np.nan
            y1_r       = np.nan
        else:
            radius_inv = 1. / radius
            if lame_2shear == 0.0:
                dy1_dr = np.nan
            else:
                dy1_dr = (1. / lame_2shear) * (y2 - (lame * radius_inv) * (2. * y1 - degree_l * (degree_l + 1.) * y3))
            y3_r       = y3 * radius_inv
            y1_r       = y1 * radius_inv

        for ci in range(n_colatitude):
            colatitude = colatitude_array[ci]
            sin_theta = np.sin(colatitude)
            if sin_theta == 0.0:
                sin_theta_inv = np.nan
            else:
                sin_theta_inv = 1. / sin_theta
            tan_theta = np.tan(colatitude)
            if tan_theta == 0.0:
                cot_theta = np.nan
            else:
                cot_theta = 1. / tan_theta
            for li in range(n_longitude):
                # longitude = longitude_domain[li]
                for ti in range(n_time):
                    # time = time_domain[ti]
                    tp       = tidal_potential[li, ci, ti]
                    tp_p_t   = tidal_potential_partial_theta[li, ci, ti]
                    tp_p_p   = tidal_potential_partial_phi[li, ci, ti]
                    tp_p2_tp = tidal_potential_partial2_theta_phi[li, ci, ti]
                    tp_p2_t2 = tidal_potential_partial2_theta2[li, ci, ti]
                    tp_p2_p2 = tidal_potential_partial2_phi2[li, ci, ti]

                    # Any optimizations
                    y1_r_tidal_potential = y1_r * tp
                    s2_t1 = (sin_theta_inv**2) * tp_p2_p2 + cot_theta * tp_p_t
                    s4_t0 = tp_p_p * sin_theta_inv
                    s5_t0 = 2. * (tp_p2_tp - cot_theta * tp_p_p) * sin_theta_inv

                    # \epsilon_{rr}
                    strains[0, ri, li, ci, ti] = dy1_dr * tp
                    # \epsilon_{\theta\theta}
                    strains[1, ri, li, ci, ti] = y3_r * tp_p2_t2 + y1_r_tidal_potential
                    # \epsilon_{\phi\phi}
                    # There is a typo in Tobie+2005 with in both the \theta\phi and \phi\phi components of the strain tensor,
                    #    as pointed out in Kervazo et al (2021; A&A) Appendix D
                    strains[2, ri, li, ci, ti] = y1_r_tidal_potential + y3_r * s2_t1

                    # TODO: "Note that the non-diagonal strains (eh/, e/r, erh) in Takeuchi and Saito (1972)
                    #  must be multiplied by 1/2 if they are to be the components of the strain tensor" Is this right?
                    #  adding in the 1/2s now...
                    # \epsilon_{r\theta}
                    strains[3, ri, li, ci, ti] = y4_shear * tp_p_t / 2.
                    # \epsilon_{r\phi}
                    strains[4, ri, li, ci, ti] = y4_shear * s4_t0 / 2.
                    # \epsilon_{\theta\phi}
                    strains[5, ri, li, ci, ti] = y3_r * s5_t0 / 2.

                    # Calculate stress assuming isotropic media (Takeuchi & Saito (1972))
                    # First calculate the trace of strain \epsilon_{rr} + \epsilon_{\theta\theta} + \epsilon_{\phi\phi}
                    strain_trace = strains[0, ri, li, ci, ti] + strains[1, ri, li, ci, ti] + strains[2, ri, li, ci, ti]
                    strain_trace_lame = lame * strain_trace

                    # Constitutive equation: \sigma_{ij} = 2\bar{\mu} \epsilon_{ij} + \bar{\lambda}\epsilon_{kk}\delta_{ij}
                    for k in range(6):
                        if k < 3:
                            stresses[k, ri, li, ci, ti] = (2. * shear) * strains[k, ri, li, ci, ti] + strain_trace_lame
                        else:
                            stresses[k, ri, li, ci, ti] = (2. * shear) * strains[k, ri, li, ci, ti]

    return strains, stresses
