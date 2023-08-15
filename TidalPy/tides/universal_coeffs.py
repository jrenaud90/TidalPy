""" Universal coefficients used to calculate tidal heating and tidal potential derivatives
    Precomputed here to avoid calls to gamma functions in expensive loops.
    Defined as:
        [(l - m)!  / (l + m)!] * (2 - delta_0m)
    Stored as:
        [order_l] [m]
        Starting with order_l = 2 and m = 0
"""

from TidalPy.exceptions import TidalPyValueException
from TidalPy.utilities.performance.numba import njit


@njit(cacheable=True)
def get_universal_coeffs(order_l: int):
    # TODO: Right now this is defined inside the function to ensure that it is compiled correctly by njit - if we make a typed dict it may be possible to pull it outside the function for better optimization
    universal_coeffs_by_orderl_minus2 = (
        # l = 2
        {
            0: 1.,
            1: 1. / 3.,
            2: 1. / 12.
            },
        # l = 3
        {
            0: 1.,
            1: 1. / 6.,
            2: 1. / 60.,
            3: 1. / 360.
            },
        # l = 4
        {
            0: 1.,
            1: 1. / 10.,
            2: 1. / 180.,
            3: 1. / 2520.,
            4: 1. / 20160.
            },
        # l = 5
        {
            0: 1.,
            1: 1. / 15.,
            2: 1. / 420.,
            3: 1. / 10080.,
            4: 1. / 181440.,
            5: 1. / 1814400.
            },
        # l = 6
        {
            0: 1.,
            1: 1. / 21.,
            2: 1. / 840.,
            3: 1. / 30240.,
            4: 1. / 907200.,
            5: 1. / 19958400.,
            6: 1. / 239500800.
            },
        # l = 7
        {
            0: 1.,
            1: 1. / 28.,
            2: 1. / 1512.,
            3: 1. / 75600,
            4: 1. / 3326400.,
            5: 1. / 119750400.,
            6: 1. / 3113510400.,
            7: 1. / 43589145600.
            }
        )

    if order_l < 2:
        raise TidalPyValueException('Tidal order l must be an integer >= 2.')

    return universal_coeffs_by_orderl_minus2[order_l - 2]
