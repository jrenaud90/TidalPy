""" Universal coefficients used to calculate tidal heating and tidal potential derivatives

    Precomputed here to avoid calls to gamma functions in expensive loops.

    Defined as:
        [(l - m)!  / (l + m)!] * delta(2 - delta_0m)

    Stored as:
        [order_l] [m]

"""

universal_coeffs_byorderl = (
    (
        1.,
        2. / 6.,
        2. / 24.
    ),
    (
        1.,
        4. / 24.,
        2. / 120.,
        2. / 720.
    ),
    (
        1.,
        12. / 120.,
        4. / 720.,
        2. / 5040.,
        2. / 40320.
    ),
    (
        1.,
        48. / 720.,
        12. / 5040.,
        4. / 40320.,
        2. / 362880.,
        2. / 3628800.
    )
)