import numpy as np

tide_defaults = {
    'base': {
        'eccentricity_truncation_lvl':     6,
        'max_tidal_order_l':               2,
        'obliquity_tides_on':              True,
        'use_planet_params_for_love_calc': True
    },
    'simple': {
        'fixed_q':                         9000.,
        'static_k2':                       .37,
        # fixed_dt_coeff used for CTL method
        'fixed_dt_coeff':                  (100. * (2. * np.pi / (86400. * 1.)))**(-1),
        'eccentricity_truncation_lvl':     6,
        'max_tidal_order_l':               2,
        'obliquity_tides_on':              True,
        'use_planet_params_for_love_calc': True,
        'use_ctl':                         False
    },
    'layered': {
        'eccentricity_truncation_lvl':     6,
        'max_tidal_order_l':               2,
        'obliquity_tides_on':              True,
        'use_planet_params_for_love_calc': False
    }
}
