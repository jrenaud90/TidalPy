## OLD G FUNC CALCULATIONS
# Newcomb Operators
@lru_cache(maxsize=200)
def newcomb(a, b, c, d):
    if c < 0 or d < 0:
        return 0
    if d > c:
        return newcomb(a=a, b=-b, c=d, d=c)
    if c == 0 and d == 0:
        return 1
    if c == 1 and d == 0:
        return b - Rational(a, 2)
    if c == 0 and d == 1:
        return -b - Rational(a, 2)
    if d == 0:
        nc = Rational(1, (4 * c)) * (2 * (2 * b - a) * newcomb(a, b + 1, c - 1, 0) +
                                     (b - a) * newcomb(a, b + 2, c - 2, 0))
        return nc
    if d != 0:
        summation_term = 0
        j = 2
        while True:
            if (c - j) < 0 or (d - j) < 0:
                break
            next_newcomb = newcomb(a, b, c - j, d - j)
            coeff = binomial_coeff(1.5, j)
            summation_term += coeff * ((-1)**j) * next_newcomb
            j += 1

        nc = Rational(1, (4 * d)) * ((-2 * (2 * b + a)) * newcomb(a, b - 1, c, d - 1) -
                                     (b + a) * newcomb(a, b - 2, c, d - 2) -
                                     (c - 5 * d + 4 + 4 * b + a) * newcomb(a, b, c - 1, d - 1) +
                                     (2 * (c - d + b)) * summation_term)
        return nc

    raise ValueError('How did you get here?')


# Hansen Coefficients
@lru_cache(maxsize=200)
def hansen(a, b, c, eccentricity, cutoff_power):
    outer_power = abs(c - b)
    alpha = max(0, c - b)
    beta = max(0, b - c)
    summation = 0
    sigma = 0
    while True:
        expo = int(2 * sigma)
        if (expo + outer_power) > cutoff_power:
            break

        newcomb_val = newcomb(a, b, sigma + alpha, sigma + beta)
        print(newcomb_val)
        summation += eccentricity**(expo) * newcomb_val
        sigma += 1

    return (eccentricity**outer_power) * summation


@lru_cache(maxsize=200)
def hansen_2(a, b, c, eccentricity, cutoff_power):
    cutoff_power_touse = cutoff_power + 4

    beta = eccentricity / (1 + (1 - eccentricity**2)**(1 / 2))
    beta = taylor(beta, eccentricity, cutoff_power_touse)

    large_num = 1000000
    max_s = a - b + 1
    max_t = a + b + 1
    if max_s < 0:
        max_s = large_num
    if max_t < 0:
        max_t = large_num
    print(max_s, max_t)

    summation = 0
    for s in range(0, max_s + 1):
        if s == large_num:
            raise ValueError('Large Number (s-sum) Exceeded.')
        for t in range(0, max_t + 1):
            if t == large_num:
                raise ValueError('Large Number (t-sum) Exceeded.')
            if (s + t > cutoff_power_touse):
                break
            coeff_1 = binomial_coeff(a - b + 1, s)
            coeff_2 = binomial_coeff(a + b + 1, t)
            e_term = (-beta)**(s + t)
            bess = besselj_func(c - b - s + t, c * eccentricity, cutoff_power_touse)
            summation += coeff_1 * coeff_2 * e_term * bess
        else:
            continue
        break

    res = (1 + beta**2)**(-a - 1) * summation
    return taylor(res, eccentricity, cutoff_power)


# Eccentricity Functions
@lru_cache(maxsize=50)
def G_func_1(l, p, q, eccentricity, cutoff_power, hfunc=1):
    """ Defined by Murray & Dermott """
    if hfunc == 1:
        return hansen(a=-l - 1, b=l - 2 * p, c=l - 2 * p + q, eccentricity=eccentricity, cutoff_power=cutoff_power)
    elif hfunc == 2:
        return hansen_2(a=-l - 1, b=l - 2 * p, c=l - 2 * p + q, eccentricity=eccentricity, cutoff_power=cutoff_power)
    else:
        return hansen_3(a=-l - 1, b=l - 2 * p, c=l - 2 * p + q, eccentricity=eccentricity, cutoff_power=cutoff_power)


@lru_cache(maxsize=50)
def G_func_2(l, p, q, eccentricity, cutoff_power):
    """ See appendix in Veras et al. 2019. Defined using Bessel functions. """

    cutoff_power_touse = cutoff_power + 1
    #     max_s = 2*p + 1
    #     max_u = 2*l - 2*p + 1
    if p - l >= 0:
        max_s = 2 * p - 2 * l
    else:
        max_s = 1000000
    if p <= 0:
        max_u = -2 * p
    else:
        max_u = 10000000000

    beta = eccentricity / (1 + (1 - eccentricity**2)**(1 / 2))
    beta = taylor(beta, eccentricity, cutoff_power_touse + 2)
    print(max_s, max_u)
    summation = 0
    for s in range(0, max_s + 1):
        for u in range(0, max_u + 1):
            if (s + u > cutoff_power_touse + 2):
                break
            coeff_1 = binomial_coeff(2 * p - 2 * l, s)
            coeff_2 = binomial_coeff(-2 * p, u)
            e_term = (-beta)**(s + u)
            bess = besselj_func(q - s + u, (l - 2 * p + q) * eccentricity, cutoff_power_touse)
            summation += coeff_1 * coeff_2 * e_term * bess
        else:
            continue
        break

    res = summation * (1 + beta**2)**(l)
    return taylor(res, eccentricity, cutoff_power_touse).removeO()


if use_md_gfunc:
    G_func = G_func_1
else:
    G_func = G_func_2