from .ctl_funcs import linear_dt, linear_dt_with_q

known_ctl_methods = {
    'linear_simple'       : linear_dt,
    'linear_simple_with_q': linear_dt_with_q
    }

ctl_method_input_getters = {
    'linear_simple'       : (('tides', 'fixed_dt'),),
    'linear_simple_with_q': (('tides', 'fixed_dt'), ('tides', 'fixed_q'))
    }
