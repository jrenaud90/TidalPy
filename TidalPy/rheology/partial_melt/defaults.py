partial_melt_defaults = {
    'ice' : {
        'model'       : 'off',
        'solidus'     : 270.,
        'liquidus'    : 273.15,
        'liquid_shear': 1.e-5
    },
    'rock': {
        'model'                 : 'henning',
        'solidus'               : 1600.,
        'liquidus'              : 2000.,
        'liquid_shear'          : 1.0e-5,
        'fs_visc_power_slope'   : 27000.0,
        'fs_visc_power_phase'   : 1.0,
        'fs_shear_power_slope'  : 82000.0,
        'fs_shear_power_phase'  : 40.6,
        'crit_melt_frac'        : 0.5,
        'crit_melt_frac_width'  : 0.05,
        'hn_visc_slope_1'       : 13.5,
        'hn_visc_slope_2'       : 370.0,
        'hn_shear_param_1'      : 40000.0,
        'hn_shear_param_2'      : 25.0,
        'hn_shear_falloff_slope': 700.0
    },
    'iron': {
        'model'       : 'off',
        'solidus'     : 4000.,
        'liquidus'    : 5000.,
        'liquid_shear': 1.e-5
    }
}
