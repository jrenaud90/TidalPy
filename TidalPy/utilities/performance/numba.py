import os
from functools import wraps

import numpy as np
from scipy.special import gamma

from ... import config

use_numba_cfg = config['use_numba']
use_numba_env = True
if 'NUMBA_DISABLE_JIT' in os.environ:
    use_numba_env = (os.environ['NUMBA_DISABLE_JIT'] == 0)
use_numba = use_numba_cfg and use_numba_env
cache_numba = config['cache_numba']

if use_numba:
    import numba

    def njit(*args, **kwargs):
        if 'cacheable' in kwargs:
            if kwargs['cacheable'] and cache_numba:
                kwargs['cache'] = True
            else:
                kwargs['cache'] = False
            del kwargs['cacheable']
            return numba.njit(*args, **kwargs)
        return numba.njit(*args, **kwargs)

    vectorize = numba.vectorize
    float64 = numba.float64
    int64 = numba.int64

else:
    vectorize = np.vectorize
    float64 = np.float64
    int64 = np.int64

    def njit(*args, **kwargs):
        def njit_inner(func):

            return func
        return njit_inner

    def nbTuple(*args):
        return None

# Special functions that work with njit
# njit does not support wrapping around scipy's gamma function. Instead we use a lookup table technique. This actually
#     ends up being faster than the gamma function, but it will give incorrect values if you are:
#     1). Outside its defined range
#     2). Taking a factorial of a number that is 0.01 / 2. away from a predefined number.
#     To counter this we add a list of commonly used values of the Andrade alpha parameter, as it is usually what we are
#     taking the factorial of.
_predefined_range = np.linspace(0., 1., 1000)
_common_alphas = np.asarray([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                             0.9, 0.95])
_predefined_inputs = np.sort(np.concatenate((_predefined_range, _common_alphas)))
_factorials = gamma(_predefined_inputs + 1.0)


@njit()
def find_factorial(number: float) -> float:
    """ Find's the factorial of a number based on a look-up table. 'number' must be between 0 and 1.

    The lookup table approach allows this function to be njit'd and used in other njit functions.

    Parameters
    ----------
    number : float
        Number to be factorialized. Must be between 0 and 1

    Returns
    -------
    factorial : float
        Factorial of 'number'
    """

    index = np.abs(_predefined_inputs - number).argmin()
    factorial = _factorials[index]

    return factorial
