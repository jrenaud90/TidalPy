world_defaults = {
    'base'    : {
        'name'                        : 'unknown_world_basetype',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',  # Options are no_eccentricity, williams, or mendez
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'use_real_moi'                : True
    },
    'simple_tide'   : {
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
        'use_real_moi'                : True
    },
    'gasgiant': {
        'name'                        : 'unknown_world_gasgianttype',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True
    },
    'star'    : {
        'name'                        : 'unknown_world_startype',
        'store_tides_config_in_world' : True,
        # Most of the stuff for the star class is not actually used, but having this here allows it to share the same
        # base class as the other planet types. Perhaps one day it would be a good idea to split the base classes of
        # stars and planets so that this stuff can be left off. For now though it really doesn't hurt anything.
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True
    },
    'tidal'   : {
        'name'                        : 'unknown_world_tidaltype',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True
    },
    'layered'   : {
        'name'                        : 'unknown_world_thermaltype',
        'store_tides_config_in_world' : True,
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'fixed_time_lag'              : 1.,
        'orbital_truncation_level'    : 2,
        'tidal_order_l'               : 2,
        'use_real_moi'                : True
    }
}
