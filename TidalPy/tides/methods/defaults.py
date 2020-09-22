import numpy as np

tide_defaults = {
    'base'         : {
        'eccentricity_truncation_lvl'    : 6,
        'max_tidal_order_l'              : 2,
        'obliquity_tides_on'             : True,
        'use_planet_params_for_love_calc': True,
        'multiply_modes_by_sign'         : True,
        'slices'                         : 100,
    },
    'global_approx': {
        'fixed_q'                        : 1000.,
        'static_k2'                      : .37,
        'eccentricity_truncation_lvl'    : 6,
        'max_tidal_order_l'              : 2,
        'obliquity_tides_on'             : True,
        'use_planet_params_for_love_calc': True,
        'multiply_modes_by_sign'         : True,
        'use_ctl'                        : False,
        # ctl_calc_method and fixed_dt used for CTL method
        'ctl_calc_method'                : 'linear_simple',
        'fixed_dt'                       :  (1. / 100.) * (2. * np.pi / (86400. * 10.))**(-1),
        'slices'                         : 100,
    },
    'layered'      : {
        'eccentricity_truncation_lvl'    : 6,
        'max_tidal_order_l'              : 2,
        'obliquity_tides_on'             : True,
        'use_planet_params_for_love_calc': False,
        'multiply_modes_by_sign'         : True
    }
}
