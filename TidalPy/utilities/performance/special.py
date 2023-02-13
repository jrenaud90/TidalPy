""" Special math functions are defined here.
"""

from scipy import special
from numba import njit


@njit(cache=True)
def find_factorial(number: float) -> float:
    """ Find's the factorial of a number based on SciPy's Gamma function.

    Notes
    -----
    This uses the numba-scipy package to register overloads for scipy's special functions. This allows Numba to njit
    these functions, including gamma.

    Parameters
    ----------
    number : float
        Number to be factorialized.

    Returns
    -------
    factorial : float
        Factorial of 'number'
    """

    return special.gamma(number + 1.)
