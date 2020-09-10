world_defaults = {
    'base'    : {
        'name'                        : 'unknown_world_base_type',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',  # Options are no_eccentricity, williams, or mendez
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'use_real_moi'                : True
    },
    'simple_tidal'   : {
        'name'                        : 'unknown_world_simple_tide_type',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True,
        'tides_on'                    : True
    },
    'gas_giant': {
        'name'                        : 'unknown_world_gas_giant_type',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True,
        'tides_on'                    : True
    },
    'star'    : {
        'name'                        : 'unknown_world_star_type',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True,
        'tides_on'                    : False
    },
    'layered'   : {
        'name'                        : 'unknown_world_thermal_type',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True,
        'tides_on'                    : False
    }
}
