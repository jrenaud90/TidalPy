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


@njit(cacheable=True)
def calculate_displacements(tidal_potential: np.ndarray,
                            tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
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


# @njit(cacheable=True)
def calculate_strain(tidal_potential: np.ndarray,
                     tidal_potential_partial_theta: np.ndarray, tidal_potential_partial_phi: np.ndarray,
                     tidal_potential_partial2_theta2: np.ndarray, tidal_potential_partial2_phi2: np.ndarray,
                     tidal_potential_partial2_theta_phi: np.ndarray,
                     tidal_solution_y: np.ndarray, tidal_solution_y_derivative: np.ndarray,
                     colatitude: FloatArray,
                     radius: np.ndarray, shear_moduli: np.ndarray):
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

    Returns
    -------
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
    """

    # Optimizations
    sin_theta = np.sin(colatitude)
    cos_theta = np.cos(colatitude)
    cot_theta = cos_theta / sin_theta

    # Shortcuts
    # Method 1 - Definition:
    dy1_dr = tidal_solution_y_derivative[0]
    # Method 2 - Euler method
    # if ri == 0:
    #     dy1_dr = tidal_solution_y[0, 1:] / radius
    # else:
    #     dy1_dr = (tidal_solution_y[0, 1:] - tidal_solution_y[0, :]) / (radius - radius[ri-1])

    y1 = tidal_solution_y[0]
    y3 = tidal_solution_y[2]
    y4 = tidal_solution_y[3]

    # TB05 Eq. 10; B13 Eq. 8 & 9
    # TODO: Missing the sum over l and m used in TB05
    e_rr = dy1_dr * tidal_potential
    e_rth = y4 * tidal_potential_partial_theta / shear_moduli
    e_rph = y4 * tidal_potential_partial_phi / (shear_moduli * sin_theta)

    e_thth = (1. / radius) * (y3 * tidal_potential_partial2_theta2 + y1 * tidal_potential)
    e_phph = (1. / radius) * \
             (y1 * tidal_potential + (y3 / sin_theta**2) *
              (tidal_potential_partial2_phi2 + cos_theta * sin_theta * tidal_potential_partial_theta))
    e_thph = (2. * y3 / (radius * sin_theta)) * \
              (tidal_potential_partial2_theta_phi - cot_theta * tidal_potential_partial_phi)

    # The (1/2) in the off-diagonal terms are due to these components appearing twice in the tensor
    e_rth *= (1. / 2.)
    e_rph *= (1. / 2.)
    e_thph *= (1. / 2.)

    # Build strain tensor matrix. Off-diagonal terms are divided by 2 in the strain tensor.
    # # TODO: Currently njit does not like the np.asarray() - perhaps in the future... (Remove 1/2 above)
    # strain_tensor = np.asarray([[e_rr, e_rth / 2., e_rph / 2.],
    #                             [e_rth / 2., e_thth, e_thph / 2.],
    #                             [e_rph / 2., e_thph / 2., e_phph]], dtype=np.complex128)

    return e_rr, e_thth, e_phph, e_rth, e_rph, e_thph