import math
from functools import lru_cache

from scipy.special import gamma
from sympy import Rational

from .general_math import binomial_coeff, besselj_func
from .sympy_help import taylor

from time import time

@lru_cache(maxsize=10)
def calc_bessel_beta(eccentricity, cutoff_power):
    beta = (1 - (1 - eccentricity**2)**(1 / 2)) / eccentricity
    beta = taylor(beta, eccentricity, cutoff_power).removeO()
    return beta

@lru_cache(maxsize=500)
def hansen_bessel(a, b, c, eccentricity, cutoff_power):

    beta = calc_bessel_beta(eccentricity, cutoff_power)

    progress = 0
    outer_sum = 0
    p = 0
    terms = 0
    while True:
        if p > cutoff_power + 1:
            break

        inner_sum = 0
        for h in range(0, p + 1):
            coeff_1 = binomial_coeff(a + b + 1, p - h)
            coeff_2 = binomial_coeff(a - b + 1, h)
            bess = besselj_func(c - b + p - 2 * h, c * eccentricity, cutoff_power)
            inner_sum += coeff_1 * coeff_2 * bess
            terms += 1

        outer_sum += inner_sum * (-beta)**p

        # If all of the terms are simply added then sympy has a real hard time taylor expanding at the end
        #    (very heavy computation time as cutoff_power > 10). To help with this, let us pick some interval where we
        #    stop the calculation and perform an early taylor expansion so that the number of terms does not grow so
        #    large for the final taylor series.
        if terms >= 20:
            progress += outer_sum
            progress = taylor(progress, eccentricity, cutoff_power)
            outer_sum = 0
            terms = 0

        p += 1

    # Depending on the number of terms, there may be a bit left over in the outer sum - this will grab that.
    if outer_sum != 0:
        progress += outer_sum

    res = (1 + beta**2)**(-a - 1) * progress
    res = taylor(res, eccentricity, cutoff_power)

    return res

@lru_cache(maxsize=100)
def hansen_kIsZero_nIsPos(n, m, eccentricity, cutoff_power, force_break=True):
    """

    k = 0, n >= 0

    See: Laskar & Boue (2010)

    """

    if n < 0:
        raise ValueError('Incorrect Domain')

    is_exact = True
    outer_coeff = (-1**m) * Rational(gamma(1 + n + m + 1), gamma(1 + n + 1))

    summation = 0

    # TODO: is the floor correct here?
    top_sum = math.floor((1 + n - m) / 2)
    for j in range(0, top_sum + 1):
        if force_break and (m + 2 * j) > cutoff_power:
            is_exact = False
            break
        coeff = Rational(gamma(1 + n - m + 1), gamma(j + 1) * gamma(m + j + 1) * gamma(1 + n - m - 2 * j + 1))
        eccen = (eccentricity / 2)**(m + 2 * j)
        summation += coeff * eccen
    return is_exact, outer_coeff * summation

@lru_cache(maxsize=100)
def hansen_kIsZero_nIsNeg(n, m, eccentricity, cutoff_power, force_break=True):
    """

    k = 0, n < -1

    See: Laskar & Boue (2010)

    """

    if n >= -1:
        raise ValueError('Incorrect Domain')

    if m < 0:
        # For k = 0, Hansen(n, m) = Hansen(n, -m)
        return hansen_kIsZero_nIsNeg(n, -m, eccentricity, cutoff_power, force_break=force_break)

    is_exact = True
    n_pos = -n
    if m >= n_pos - 1:
        # Hansen Coefficients have the property that they are equal to 0 for m > |n| (Hughes 1981); see also L&B 2010
        return is_exact, 0
    else:
        outer_coeff = ((1 - eccentricity**2)**(Rational(3, 2) - n_pos))

        summation = 0
        # TODO: is the floor correct here?
        top_sum = math.floor((n_pos - 2 - m) / 2)
        for j in range(0, top_sum + 1):
            if force_break and (m + 2 * j) > cutoff_power:
                is_exact = False
                break
            coeff = Rational(gamma(n_pos - 2 + 1), gamma(j + 1) * gamma(m + j + 1) * gamma(n_pos - 2 - m - 2 * j + 1))
            eccen = (eccentricity / 2)**(m + 2 * j)
            summation += coeff * eccen
        return is_exact, outer_coeff * summation

@lru_cache(maxsize=200)
def hansen_wrapper(n, m, k, eccentricity, cutoff_power,
                   force_break: bool = True):
    # Some hansen numbers can be provided with exact precision (if k==0)
    is_exact = False

    if k == 0:
        is_exact = True
        if m < 0:
            return hansen_wrapper(n, -m, 0, eccentricity, cutoff_power,
                                  force_break=force_break)

        # When k==0 then the Hansen coefficients are exact (see Laskar and Boue 2010)
        if n == -1:
            # Case where n == -1
            if m == 0:
                return is_exact, 1
            elif m == 1:
                return is_exact, ((1 - eccentricity**2)**(1 / 2) - 1) / eccentricity
            else:
                return is_exact, 0
        elif n < -1:
            # Case where n < -1
            return hansen_kIsZero_nIsNeg(n, m, eccentricity, cutoff_power, force_break=force_break)

        else:
            # Case where n >= 0
            return hansen_kIsZero_nIsPos(n, m, eccentricity, cutoff_power, force_break=force_break)

    else:
        # k != 0 is not exact and requires a truncation on a series (see Renaud et al. 2020)
        if m <= 0 and k < 0:
            # TODO: Where is this from?
            return is_exact, hansen_bessel(n, -m, -k, eccentricity, cutoff_power)

        return is_exact, hansen_bessel(n, m, k, eccentricity, cutoff_power)
