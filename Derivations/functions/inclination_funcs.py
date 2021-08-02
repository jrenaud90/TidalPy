from functools import lru_cache

from scipy.special import gamma
import math
from sympy import Symbol, Rational, expand_trig, simplify, trigsimp
from sympy.functions.elementary.trigonometric import sin, cos

from .general_math import binomial_coeff
from .sympy_help import taylor

sin_f_S = Symbol('S___S', real=True)
sin_f_H = Symbol('S___H', real=True)
cos_f_S = Symbol('C___S', real=True)
cos_f_H = Symbol('C___H', real=True)


## Inclination Functions
@lru_cache(maxsize=50)
def F_func(l, m, p, inclination, cut_off_power=None, auto_taylor: bool = False, run_trigsimp: bool = True):
    lower_sum_bound = max(0, l - m - 2 * p)
    upper_sum_bound = min(l - m, 2 * l - 2 * p)

    cos_f = cos(inclination / 2)
    sin_f = sin(inclination / 2)

    # See the discussion in Gooding & Wagner 2008 and in the appendix of Renaud+2020. also Veras et al 2019
    outer_coeff = Rational(gamma(l + m + 1), (2**l * gamma(p + 1) * gamma(l - p + 1)))

    summation = 0.
    for lam in range(lower_sum_bound, upper_sum_bound + 1):
        cos_expo = int(3 * l - m - 2 * p - 2 * lam)
        sin_expo = int(m - l + 2 * p + 2 * lam)

        term_1 = (-1)**lam
        term_2 = binomial_coeff(2 * l - 2 * p, lam)
        term_3 = binomial_coeff(2 * p, l - m - lam)
        term_4 = cos_f**cos_expo
        term_5 = sin_f**sin_expo

        summation += term_1 * term_2 * term_3 * term_4 * term_5

    result = outer_coeff * summation

    # It appears that this entier function should be used. This is a problem when F is used in place of F^2
    # See discussion around Eq. 2 of Gooding and Wagner 2010
    int_func = (1. / 2.) * (l - m + 1.)
    # Using the Entier functions "Integer part" defined as floor if x >= 0 and ceiling otherise
    if int_func >= 0:
        ent_scale = math.floor(int_func)
    else:
        ent_scale = math.ceil(int_func)
    result = (-1)**ent_scale * result
    if auto_taylor:
        if cut_off_power is None:
            raise Exception('Cutoff Power Required')
        result = taylor(result, inclination, cut_off_power + 1).removeO()
    elif run_trigsimp:
        result = simplify(expand_trig(trigsimp(result)))

    return result
