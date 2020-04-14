""" Universal coefficients used to calculate tidal heating and tidal potential derivatives
    Precomputed here to avoid calls to gamma functions in expensive loops.
    Defined as:
        [(l - m)!  / (l + m)!] * (2 - delta_0m)
    Stored as:
        [order_l] [m]
        Starting with order_l = 2 and m = 0
"""
from ..exceptions import TidalPyValueException
from ..performance import njit

@njit
def get_universal_coeffs(order_l):

    # TODO: Right now this is defined inside the function to ensure that it is compiled correctly by njit - if we make a typed dict it may be possible to pull it outside the function for better optimization
    universal_coeffs_byorderl = (
        {
            0: 1.,
            1: 2. / 6.,
            2: 2. / 24.
        },
        {
            0: 1.,
            1: 4. / 24.,
            2: 2. / 120.,
            3: 2. / 720.
        },
        {
            0: 1.,
            1: 12. / 120.,
            2: 4. / 720.,
            3: 2. / 5040.,
            4: 2. / 40320.
        },
        {
            0: 1.,
            1: 48. / 720.,
            2: 12. / 5040.,
            3: 4. / 40320.,
            4: 2. / 362880.,
            5: 2. / 3628800.
        }
    )

    if order_l < 2:
        raise TidalPyValueException('Tidal order l must be an integer >= 2.')

    return universal_coeffs_byorderl[order_l - 2]