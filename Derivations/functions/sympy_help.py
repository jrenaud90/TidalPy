from pprint import pprint
from sympy import series, nsimplify, Symbol, init_printing, Function, collect
from sympy.functions.special.polynomials import assoc_legendre
from IPython.display import display

USE_PRETTY_PRINT = True

# Symbolic Functions and Variables
time_var = Symbol('t', positive=True, real=True)
eccentricity = Symbol('e', positive=True, real=True)
obliquity_sat = Symbol('I___S', positive=True, real=True)
obliquity_host = Symbol('I___H', positive=True, real=True)

mean_n = Symbol('n', real=True, positive=True)
spin_host = Symbol('Omega___H', real=True, positive=True)
spin_sat = Symbol('Omega___S', real=True, positive=True)
moi_host = Symbol('C___H', real=True, positive=True)
moi_sat = Symbol('C___S', real=True, positive=True)
r_sat = Symbol('r___S', real=True, positive=True)
r_host = Symbol('r___S', real=True, positive=True)

love2_num_host = Function('Xi___H2', real=True)
love2_num_sat = Function('Xi___S2', real=True)
love3_num_host = Function('Xi___H3', real=True)
love3_num_sat = Function('Xi___S3', real=True)
love4_num_host = Function('Xi___H4', real=True)
love4_num_sat = Function('Xi___S4', real=True)
love5_num_host = Function('Xi___H5', real=True)
love5_num_sat = Function('Xi___S5', real=True)
love6_num_host = Function('Xi___H6', real=True)
love6_num_sat = Function('Xi___S6', real=True)
love7_num_host = Function('Xi___H7', real=True)
love7_num_sat = Function('Xi___S7', real=True)
love8_num_host = Function('Xi___H8', real=True)
love8_num_sat = Function('Xi___S8', real=True)

legendreP = assoc_legendre
legendreP_symbolic = Function('P^{l}', real=True)

love_funcs_host = (love2_num_host, love3_num_host, love4_num_host, love5_num_host, love6_num_host, love7_num_host,
                   love8_num_host)
love_funcs_sat = (love2_num_sat, love3_num_sat, love4_num_sat, love5_num_sat, love6_num_sat, love7_num_sat,
                   love8_num_sat)

periapsis = Symbol('\\omega', real=True, positive=True)
mean_anon = Symbol('\\mathcal{M}', real=True, positive=True)
node = Symbol('\\Omega', real=True, positive=True)
sidereal_time = Symbol('\\Theta', real=True, positive=True)
longitude = Symbol('\\phi', real=True, positive=True)
latitude = Symbol('\\theta', real=True)
def find_varpi(l: int, m: int, p: int, q: int):

    return (l - 2 * p) * periapsis + (l - 2 * p + q) * mean_anon + m * node

love_func_lookup = {
    'host': {
        2: love2_num_host,
        3: love3_num_host,
        4: love4_num_host,
        5: love5_num_host,
        6: love6_num_host,
        7: love7_num_host,
        8: love8_num_host,
    },
    'sat': {
        2: love2_num_sat,
        3: love3_num_sat,
        4: love4_num_sat,
        5: love5_num_sat,
        6: love6_num_sat,
        7: love7_num_sat,
        8: love8_num_sat,
    }
}

sign_func = Function('Upsilon')

# Planetary and Orbital Variables
mass_host = Symbol('M___H', positive=True, real=True)
mass_sat = Symbol('M___S', positive=True, real=True)
radius_host = Symbol('R___H', positive=True, real=True)
radius_sat = Symbol('R___S', positive=True, real=True)
semi_major_axis = Symbol('a', positive=True, real=True)
newton_G = Symbol('G', positive=True, real=True)

# Script Adjustments
if USE_PRETTY_PRINT:
    init_printing()
    disp = lambda text: display(text)
else:
    disp = print

# Helper Functions
def taylor(func, x, n, x0=0., simplify=True):
    res = series(func, x, x0, n)
    if simplify:
        res = nsimplify(res)
    return res


## Expression Functions
def sync_spin(expr, host: bool, sat: bool, max_q):

    if host:
        # Set the spin equal to the mean motion
        expr = expr.subs(spin_host, mean_n)

    if sat:
        expr = expr.subs(spin_sat, mean_n)

    # Clean up the function
    expr = clean_up(expr, max_q)

    return expr

def make_orbit_boring(expr, reduce_e):
    if reduce_e:
        expr = expr.subs(eccentricity, 0)
    expr = expr.subs(obliquity_sat, 0)
    expr = expr.subs(obliquity_host, 0)
    return expr.removeO()

def clean_up(expr, max_q):

    for i in range(max_q):
        # Pull out negatives from the love functions
        for love_func in love_funcs_host:
            expr = expr.subs(love_func(-i * mean_n), -love_func(i * mean_n))
            expr = expr.subs(love_func(-i * spin_host), -love_func(i * spin_host))

        for love_func in love_funcs_sat:
            expr = expr.subs(love_func(-i * mean_n), -love_func(i * mean_n))
            expr = expr.subs(love_func(-i * spin_sat), -love_func(i * spin_sat))

    # Set love numbers at 0 frequency to 0
    for love_func in love_funcs_host:
        expr = expr.subs(love_func(0), 0)

    for love_func in love_funcs_sat:
        expr = expr.subs(love_func(0), 0)

    expr = expr.subs(sign_func(0), 0)
    expr = expr.subs(sign_func(-mean_n), -1 * sign_func(mean_n))

    # Remove any Orders, if present
    expr = expr.removeO()
    return expr

def disp_wclean(expr, max_q):
    disp(clean_up(expr, max_q))

def flip_a2n(expr):
    return expr.subs(semi_major_axis, (newton_G * (mass_host + mass_sat) / mean_n**2)**(1 / 3))

def remove_G(expr):
    return expr.subs(newton_G, semi_major_axis**3 * mean_n**2 / (mass_host + mass_sat))

def collect_ei(expr):
    expr = collect(expr, eccentricity)
    expr = collect(expr, obliquity_host)
    expr = collect(expr, obliquity_sat)

    return expr
