""" Functions to decompose the tidal solution into useful properties.

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005), DOI: 10.1016/j.icarus.2005.04.006
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
"""

import numpy as np

from ...utilities.performance import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def decompose(tidal_y: np.ndarray, tidal_y_derivative: np.ndarray, radius_array: np.ndarray, gravity_array: np.ndarray,
              complex_shear_modulus: np.ndarray, bulk_modulus: FloatArray, order_l: int = 2):
    """ Decomposes the tidal solution (y) into useful properties.

    Note: Both the radius and tidal_y must be provided with length N (number of shells), but shear and bulk modulus
        must be provided a N-1 (skipping the innermost index)

    References
    ----------
    TB05: Eq. 33

    Parameters
    ----------
    tidal_y : np.ndarray
        Tidal propagation solution (Matrix: 6 x N) found via the propagation technique.
    tidal_y_derivative : np.ndarray
        Derivatives of the tidal propagation solutions with respect to radius (Matrix: 6 x N)
    radius_array : np.ndarray
        Radii of the world (N) [m]
    gravity_array : np.ndarray
        Acceleration due to gravity at the top of each shell within the world [m s-2].
    complex_shear_modulus : np.ndarray
        The complex shear modulus as found by the world or layer's rheology (N-1) [Pa].
    bulk_modulus : FloatArray
        The bulk modulus of the world [Pa]. Can be provided as a float (constant for the whole world) or as an array.
        If provided as an array, then it must have the same shape as the radius array (N-1).
        Note from ID:
            "K=200.0e9; // Bulk modulus (g cm-2 s-2), arbitrarily higher than K_rock (39-133 GPa) and K_ice (10.7 GPa)
            for consistency with incompressible prop mtx."
    order_l : int = 2
        Tidal Harmonic Degree Index

    Returns
    -------
    radial_sensitivity_to_shear : np.ndarray
        Radial sensistivity to shear stress as defined in TB05, Eq. 33.
    love_numbers : Tuple[np.ndarray, np.ndarray, np.ndarray]
        Tuple of the world's Love and Shida numbers: (k_numbers, h_numbers, l_numbers).

    """

    # Shortcuts
    radius = radius_array
    y1 = tidal_y[0, :]
    y2 = tidal_y[1, :]
    y3 = tidal_y[2, :]
    y4 = tidal_y[3, :]
    y5 = tidal_y[4, :]

    # As was done in ID, we invert y2 and y3 here... TODO: Check with Marc why this is.
    # Optimizations
    llp1 = order_l * (order_l + 1)
    y1mllp1y3 = 2. * y1 - order_l * (order_l + 1.) * y3

    # The gradient of y1 is needed. ID used the following method to make that calculation. However,
    #    I (JPR) believe that the derivative is provided by the definition of the tidal solutions which we can easily
    #    calculate via a matrix operation in propagate.py TODO
    #    ----ID Method----
    #        Calculate the gradient of y1
    #        dy1_dr_conj = (np.conj(y1_full[1:]) - np.conj(y1_full[:-1])) / (radius_array[1:] - radius_array[:-1])
    dy1_dr_conj = np.conj(tidal_y_derivative[0, :])

    # Radial sensitivity to shear modulus (TB05 Eq. 33)
    #     ID flips y2 and y3 are inverted here, but we have already done that in propagate.py
    shear_term1 = \
        (4. / 3.) * radius**2 / (np.abs(bulk_modulus + (4. / 3.) * complex_shear_modulus)**2) * \
        (np.abs(y2 - (((bulk_modulus - (2. / 3.) * complex_shear_modulus) / radius) * y1mllp1y3) )**2)

    shear_term2 = \
        - (4. / 3.) * radius * np.real(dy1_dr_conj * y1mllp1y3)

    shear_term3 = \
        (1. / 3.) * np.abs(y1mllp1y3)**2 + \
        llp1 * radius**2 * np.abs(y4)**2 / np.abs(complex_shear_modulus)**2

    shear_term4 = \
        order_l * (order_l**2 - 1.) * (order_l + 2.) * np.abs(y3)**2

    radial_sensitivity_to_shear = shear_term1 + shear_term2 + shear_term3 + shear_term4

    # Calculate Love and Shida numbers
    # Notes from ID:
    #    Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), opposite convention of Tobie et al. (2005, eqs. 9 & 36)
    #    And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9
    ## Okay, some clarification on this. It looks like VS04 that HH14 is based on used a different convention for y5,
    #    Tobie05's y5 = -y5 of SV; we follow that format here.
    k_numbers = y5 - 1.
    h_numbers = y1 * gravity_array
    l_numbers = y3 * gravity_array
    love_numbers = (k_numbers, h_numbers, l_numbers)

    return radial_sensitivity_to_shear, love_numbers