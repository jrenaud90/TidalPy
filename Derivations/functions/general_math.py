from functools import lru_cache

from scipy.special import gamma
from sympy import Rational, I


@lru_cache(maxsize=2000)
def binomial_coeff(top, bot):
    if top == bot:
        return 1

    if bot == 0:
        return 1

    if bot < 0:
        return 0

    if top < 0:
        if bot >= 0:
            return (-1)**bot * binomial_coeff(-top + bot - 1, bot)
        elif bot <= top:
            return (-1)**(top - bot) * binomial_coeff(-bot - 1, top - bot)
        else:
            return 0
    else:
        if 0 <= bot < top:
            return Rational(gamma(top + 1), (gamma(bot + 1) * gamma(top - bot + 1)))
        else:
            return 0


@lru_cache(maxsize=5000)
def besselj_func(a, x, cutoff):
    abs_a = abs(a)
    summation = 0
    m = 0
    while True:
        expo = 2 * m + abs_a
        if expo > cutoff:
            break

        this_term = (-1)**m * Rational(1, gamma(m + 1) * gamma(m + abs_a + 1)) * Rational(1, 2**expo) * x**expo
        summation += this_term
        m += 1

    return summation * I**(abs_a - a)
