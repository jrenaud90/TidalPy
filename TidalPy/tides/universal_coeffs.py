""" Universal coefficients used to calculate tidal heating and tidal potential derivatives
    Precomputed here to avoid calls to gamma functions in expensive loops.
    Defined as:
        [(l - m)!  / (l + m)!] * delta(2 - delta_0m)
    Stored as:
        [order_l] [m]
        Starting with order_l = 2 and m = 0
"""

from ..performance import njit

@njit
def get_universal_coeffs(order_l):
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
    return universal_coeffs_byorderl[order_l]