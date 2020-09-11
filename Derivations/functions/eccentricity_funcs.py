import math
from functools import lru_cache

from .hansen import hansen_wrapper


# Eccentricity Functions
@lru_cache(maxsize=50)
def G_func(l, p, q, eccentricity, cutoff_power, going_to_square: bool = False, use_floats=False):

    cutoff_power_to_use = cutoff_power
    if going_to_square:
        # Fix the cutoff power to ensure that precision is maximum while keeping efficiency high.
        #     First we need to use +1 for the taylor since "cutoff" means "last one we WANT to keep"
        cutoff_power_to_use = max(cutoff_power + 1 - abs(q), 3)
        # An additional efficiency gain can be made by reducing the cutoff threshold if we assume the user is going to
        #    square the result. Since the minimum e power is equal to q then we never need terms that are > cutoff - q

    # Perform some calculations here to avoid them later on

    return hansen_wrapper(n=-l-1, m=l-2*p, k=l-2*p+q,
                          eccentricity=eccentricity, cutoff_power=cutoff_power_to_use, use_floats=use_floats)