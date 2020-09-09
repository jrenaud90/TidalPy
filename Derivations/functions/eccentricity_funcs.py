from .hansen import hansen_wrapper

# Eccentricity Functions
def G_func(l, p, q, eccentricity, cutoff_power, going_to_square: bool = True):
    return hansen_wrapper(n=-l-1, m=l-2*p, k=l-2*p+q,
                          eccentricity=eccentricity, cutoff_power=cutoff_power, going_to_square=going_to_square)