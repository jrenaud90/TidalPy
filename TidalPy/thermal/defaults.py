
# Stored in World > Layer > 'cooling'
cooling_param_defaults = {
    'ice': {
        'model': 'convection',
        'convection_alpha': 1.,
        'convection_beta': 1./3.,
        'critical_rayleigh': 1600.
    },
    'rock': {
        'model': 'convection',
        'convection_alpha': 1.,
        'convection_beta': 1./3.,
        'critical_rayleigh': 1100.
    },
    'iron': {
        'model': 'off'
    }
}

# Stored in World > Layer > 'partial_melt'
partial_melter_param_defaults = {
    'ice': {
        'model': 'off',
        'solidus': 270.,
        'liquidus': 273.15,
        # Liquid Shear is just something very small.
        'liquid_shear': 1.0e-5,
    },
    'rock': {
        'model': 'henning',
        'solidus': 1600.,
        'liquidus': 2000.,
        'liquid_shear': 1.0e-5,
        'fs_visc_power_slope': 27000.0,
        'fs_visc_power_phase': 1.0,
        'fs_shear_power_slope': 82000.0,
        'fs_shear_power_phase': 40.6,
        'crit_melt_frac': 0.5,
        'crit_melt_frac_width': 0.05,
        'hn_visc_slope_1': 13.5,
        'hn_visc_slope_2': 370.0,
        'hn_shear_param_1': 40000.0,
        'hn_shear_param_2': 25.0,
        'hn_shear_falloff_slope': 700.0
    },
    'iron': {
        'model': 'off',
        'solidus': 2200.0,
        'liquidus': 3000.0,
        'liquid_shear': 1.0e-5,
    }
}