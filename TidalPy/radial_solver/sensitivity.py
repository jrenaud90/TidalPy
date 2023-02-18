""" Functionality to calculate the sensitivity parameters described in Tobie et al (2005).

References
----------

TB05  : Tobie et al. (2005; DOI: 10.1016/j.icarus.2005.04.006)
"""

import numpy as np

from TidalPy.utilities.performance import njit

@njit(cacheable=True)
def sensitivity_to_shear(
        radial_solutions: np.ndarray, radius_array: np.ndarray,
        shear_modulus_array: np.ndarray, bulk_modulus_array: np.ndarray, order_l: int = 2
        ):
    """ Decomposes the tidal solution (y) into useful properties.

    Note: Both the radius and tidal_y must be provided with length N (number of shells), but shear and bulk modulus
        must be provided a N-1 (skipping the innermost index)

    References
    ----------
    TB05: Eq. 33

    Parameters
    ----------
    radial_solutions : np.ndarray
        Viscoelastic-Gravitational radial solutions found through either the matrix propagation or
        numerical integration techniques.
        These should follow the TS72 order convention and dimensions. With the exception that there are double
        the number of values to account for the imaginary portions. E.g., y1_real, y1_imag, ... (Matrix: 12 x N)
    radius_array : np.ndarray
        Radius at the top of each radial slice throughout the world (length N) [m]
    shear_modulus_array : np.ndarray
        The shear modulus defined at each radius slice (length N) [Pa].
        Can be real or complex valued.
    bulk_modulus_array : np.ndarray
        The bulk modulus defined at each radius slice (length N) [Pa].
        Can be real or complex valued.
        Note from ID:
            "K=200.0e9; // Bulk modulus (g cm-2 s-2), arbitrarily higher than
            K_rock (39-133 GPa) and K_ice (10.7 GPa) for consistency with incompressible prop mtx."
    order_l : int = 2
        Tidal Harmonic Degree

    Returns
    -------
    radial_sensitivity_to_shear : np.ndarray
        Radial sensitivity to shear stress as defined in TB05, Eq. 33. This is a real-valued float-array (length N)

    """

    # Basic information
    len_r = radius_array.shape[0]

    # Build arrays
    shear_sensitivity = np.empty(len_r, dtype=np.float64)

    # Optimizations
    llp1     = order_l * (order_l + 1.)
    ll2m1lp2 = order_l * (order_l * order_l - 1.) * (order_l + 2.)

    # Extract the required radial solution values, isolating the real and imaginary portions.
    for r_i in range(len_r):
        y1_real = radial_solutions[0, r_i]
        y1_imag = radial_solutions[1, r_i]
        y2_real = radial_solutions[2, r_i]
        y2_imag = radial_solutions[3, r_i]
        y3_real = radial_solutions[4, r_i]
        y3_imag = radial_solutions[5, r_i]
        y4_real = radial_solutions[6, r_i]
        y4_imag = radial_solutions[7, r_i]

        # Shear and bulk may be real or complex
        shear = shear_modulus_array[r_i]
        bulk  = bulk_modulus_array[r_i]

        # Combine radial solutions into complex numbers
        y2 = y2_real + 1.0j * y2_imag
        y1_conj = y1_real - 1.0j * y1_imag

        # Find the gradient of y1_conj.
        #   The radius step size may not be even throughout the planet if there are
        #   layers with more slices. So we will use a 2nd order approach based on np.gradient method.
        r = radius_array[r_i]
        if r_i == 0:
            # First edge point
            y1_conj_rp1 = radial_solutions[0, r_i + 1] - 1.0j * radial_solutions[1, r_i + 1]  # conj(y1) at r + 1
            dr0 = radius_array[r_i + 1] - r
            y1_gradient_conj = (y1_conj_rp1 - y1_conj) / dr0
        elif r_i == len_r - 1:
            # Last end point
            y1_conj_rm1 = radial_solutions[0, r_i - 1] - 1.0j * radial_solutions[1, r_i - 1]  # conj(y1) at r - 1
            dr0 = r - radius_array[r_i - 1]
            y1_gradient_conj = (y1_conj - y1_conj_rm1) / dr0
        else:
            # Midpoints
            y1_conj_rp1 = radial_solutions[0, r_i + 1] - 1.0j * radial_solutions[1, r_i + 1]  # conj(y1) at r + 1
            y1_conj_rm1 = radial_solutions[0, r_i - 1] - 1.0j * radial_solutions[1, r_i - 1]  # conj(y1) at r - 1
            dr0 = r - radius_array[r_i - 1]  # np.diff(r)[:-1]
            dr1 = radius_array[r_i + 1] - r  # np.diff(r)[0:]

            dra = -dr1 / (dr0 * (dr0 + dr1))
            drb = (dr1 - dr0) / (dr0 * dr1)
            drc = dr0 / (dr1 * (dr0 + dr1))

            y1_gradient_conj = dra * y1_conj_rm1 + drb * y1_conj + drc * y1_conj_rp1

        # Optimizations
        r2 = r * r
        y1y3_term_real = 2. * y1_real - llp1 * y3_real
        y1y3_term_imag = 2. * y1_imag - llp1 * y3_imag
        y1y3_term = y1y3_term_real + 1.0j * y1y3_term_imag
        y1y3_term_abs2 = (y1y3_term_real * y1y3_term_real) + (y1y3_term_imag * y1y3_term_imag)
        y4_abs2 = y4_real * y4_real + y4_imag * y4_imag
        y3_abs2 = y3_real * y3_real + y3_imag * y3_imag

        # Find shear sensitivity
        if r == 0.:
            # The shear sensitivity is not defined at r=0
            shear_sensitivity[r_i] = np.nan
        else:
            shear_sensitivity[r_i] = \
                ((4. / 3.) * r2 / (np.abs(bulk + (4. / 3.) * shear)**2) *
                 np.abs(y2 - ((bulk - (2. / 3.) * shear) / r) * y1y3_term)**2) + \
                (-(4. / 3.) * r * np.real(y1_gradient_conj * y1y3_term) + (1. / 3.) * y1y3_term_abs2) + \
                ((llp1 * r2 * y4_abs2 / (np.abs(shear)**2)) + (ll2m1lp2 * y3_abs2))

    return shear_sensitivity
