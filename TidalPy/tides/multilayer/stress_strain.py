""" Calculation of stress, strain, and displacement based on tidal solutions and tidal potential

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

import numpy as np

from ...utilities.types import FloatArray
from ...utilities.performance import njit

def calculate_displacements(tidal_potential: FloatArray,
                            tidal_potential_partial_theta: FloatArray, tidal_potential_partial_phi: FloatArray,
                            tidal_solution_y: np.ndarray, colatitude: FloatArray):
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

def calculate_strain(tidal_potential: FloatArray,
                     tidal_potential_partial_theta: FloatArray, tidal_potential_partial_phi: FloatArray,
                     tidal_potential_partial2_theta2: FloatArray, tidal_potential_partial2_phi2: FloatArray,
                     tidal_potential_partial2_theta_phi: FloatArray,
                     tidal_solution_y: np.ndarray, tidal_solution_y_derivative: np.ndarray,
                     colatitude: FloatArray,
                     radius_array: np.ndarray, shear_moduli: np.ndarray):
    """ Calculate tidal strain tensor using the tidal potential and its partial derivatives as well as the y-solution
    vector.

    Parameters
    ----------
    tidal_potential : FloatArray
        The tidal potential
    tidal_potential_partial_theta : FloatArray
        The partial derivative of the tidal potential with respect to the colatitude
    tidal_potential_partial_phi : FloatArray
        The partial derivative of the tidal potential with respect to the longitude
    tidal_potential_partial2_theta2 : FloatArray
        The 2nd partial derivative of the tidal potential with respect to the colatitude
    tidal_potential_partial2_phi2 : FloatArray
        The 2nd partial derivative of the tidal potential with respect to the longitude
    tidal_potential_partial2_theta_phi : FloatArray
        The 2nd partial derivative of the tidal potential with respect to the colatitude and the longitude
    tidal_solution_y : np.ndarray
        Matrix [6 x N] of tidal solutions often denoted as a lower-case "y"
    tidal_solution_y_derivative: np.ndarray
        Matrix [6 x N] of the derivative of the tidal solutions with respect to radius (dy/dr).
    colatitude : FloatArray
        Co-latitude where the calculations are conducted
    radius_array : np.ndarray
        Radius array [m]
    shear_moduli : np.ndarray
        Shear modulus as a function of radius [Pa]

    Returns
    -------
    strain_tensor : np.ndarray
        The symmetric strain tensor
    """

    # Optimizations
    sin_theta = np.sin(colatitude)
    cos_theta = np.cos(colatitude)
    cot_theta = cos_theta / sin_theta
    r_inv_array = (1. / radius_array)

    # Shortcuts
    dy1_dr = tidal_solution_y_derivative[0, :]
    y1 = tidal_solution_y[0, :]
    y3 = tidal_solution_y[2, :]
    y4 = tidal_solution_y[3, :]

    # TB05 Eq. 10; B13 Eq. 8 & 9 TODO: Missing the sum over l and m used in TB05
    strain_tensor_by_r = list()

    # TODO: Missing r dependence in the tidal potential
    for ri, radius in enumerate(radius_array):
        # Shortcuts
        dy1_dr = tidal_solution_y_derivative[0, ri]
        y1 = tidal_solution_y[0, ri]
        y3 = tidal_solution_y[2, ri]
        y4 = tidal_solution_y[3, ri]
        r_inv  = r_inv_array[ri]
        shear = shear_moduli[ri]

        e_rr = dy1_dr * tidal_potential
        e_rth = y4 * tidal_potential_partial_theta / shear
        e_rph = y4 * tidal_potential_partial_phi / (shear * sin_theta)

        e_thth = r_inv * (y3 * tidal_potential_partial2_theta2 + y1 * tidal_potential)
        e_phph = r_inv * (y1 * tidal_potential +
                          (y3 / sin_theta**2) * (tidal_potential_partial2_phi2 +
                                                 cos_theta * sin_theta * tidal_potential_partial_theta) )
        e_thph = (2. * r_inv * y3 / sin_theta) * (tidal_potential_partial2_theta_phi -
                                                  cot_theta * tidal_potential_partial_phi)

        # Off-diagonal terms are divided by 2 in the strain tensor.
        strain_tensor = [
            [
                # e_r_r
                e_rr,
                # e_r_theta
                e_rth / 2.,
                # e_r_phi
                e_rph / 2.,
            ],
            [
                # e_r_theta
                e_rth / 2.,
                # e_theta_theta
                e_thth,
                # e_theta_phi
                e_thph / 2.,
            ],
            [
                # e_r_phi
                e_rph / 2.,
                # e_theta_phi
                e_thph / 2.,
                # e_phi_phi
                e_phph
            ]
        ]

        strain_tensor_by_r.append(strain_tensor)

    strain_tensor = np.asarray(strain_tensor_by_r)

    return strain_tensor