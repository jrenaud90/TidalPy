""" This module provides several special functions that are specifically designed to work with TidalPy and its
dependencies (looking at you, Numba).

"""

from typing import TYPE_CHECKING

import numpy as np

from TidalPy.utilities.performance import njit, use_numba

if TYPE_CHECKING:
    from TidalPy.utilities.types import NumArray


def _sqrt_neg_python(z: 'NumArray', is_real: bool = False) -> 'NumArray':
    """ Square root - Allows for negative values

    Parameters
    ----------
    z : FloatArray
        Input value (domain is all positive and negative numbers)

    Returns
    -------
    z_sqrt : FloatArray
        Output value (range is all positive values and complex numbers)

    """

    if is_real:
        # First solve the square root assuming z is positive.
        z_sqrt_abs = np.sqrt(np.abs(z))

        # Now correct the negatives (this will either be a Boolean or an array of Booleans, depending upon input type)
        z_sqrt = (np.real(z) > 0.) * z_sqrt_abs + \
                 (np.real(z) < 0.) * z_sqrt_abs * 1.0j

    else:
        # This is a more "complex" process because the input could already be both negative AND complex.
        z_r = np.real(z)
        z_i = np.imag(z)
        quad = np.sqrt(z_r * z_r + z_i * z_i)

        real_part = np.sqrt((quad + z_r) / 2.)
        imag_part = np.sqrt((quad - z_r) / 2.)

        z_sqrt = real_part + \
                 (z_i != 0.) * imag_part * np.sign(z_i) * 1.0j + \
                 (z_i == 0.) * imag_part * 1.0j

    return z_sqrt


# Imaginary square roots
if use_numba:
    # TODO: Numba currently does not support wrapping np.lib.scimath.sqrt, so we have to define our own function.
    #    However, numba.njit of np.sqrt is about 10x faster than np.sqrt with floats and about 2x fast with arrays.
    #    So the above method is actually pretty efficient.
    sqrt_neg = njit(cacheable=True)(_sqrt_neg_python)
else:
    # Numpy already has a built in function to handle this. Use it instead
    def sqrt_neg(z, is_real: bool = False):
        return np.lib.scimath.sqrt(z)
