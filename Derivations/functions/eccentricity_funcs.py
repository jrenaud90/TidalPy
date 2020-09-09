from .hansen import hansen_wrapper
from .sympy_help import taylor

# Eccentricity Functions
def G_func(l, p, q, eccentricity, cutoff_power, going_to_square: bool = True):

    cutoff_power_to_use = cutoff_power
    if going_to_square:
        # Fix the cutoff power to ensure that precision is maximum while keeping efficiency high.
        #     First we need to use +1 for the taylor since "cutoff" means "last one we WANT to keep"
        cutoff_power_to_use = cutoff_power + 1 - abs(q)
        # An additional efficiency gain can be made by reducing the cutoff threshold if we assume the user is going to
        #    square the result. Since the minimum e power is equal to q then we never need terms that are > cutoff - q

    # Perform some calculations here to avoid them later on

    return hansen_wrapper(n=-l-1, m=l-2*p, k=l-2*p+q,
                          eccentricity=eccentricity, cutoff_power=cutoff_power_to_use)