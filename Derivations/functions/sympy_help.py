from pprint import pprint
from sympy import series, nsimplify, Symbol, init_printing

USE_PRETTY_PRINT = True

eccentricity = Symbol('e', positive=True, real=True)
obliquity_sat = Symbol('I__S', positive=True, real=True)
obliquity_host = Symbol('I__H', positive=True, real=True)

symbols = {
    'eccentricity': eccentricity,
    'obliquity_host': obliquity_host,
    'obliquity_sat': obliquity_sat
}

# Script Adjustments
if USE_PRETTY_PRINT:
    init_printing()
    disp = lambda text: pprint(text)
else:
    disp = print

# Helper Functions
def taylor(func, x, n, x0=0., debug=False, simplify=True):
    if debug:
        return func
    else:
        res = series(func, x, x0, n)
        if simplify:
            res = nsimplify(res)
    return res