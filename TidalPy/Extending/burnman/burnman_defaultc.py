default_burnman_configs = {
    'worlds': {
        'types': {
            'burnman': {
                'name'                                : 'unknown_world_burnman_type',
                'store_tides_config_in_world'         : True,
                'force_spin_sync'                     : True,
                'equilibrium_insolation_model'        : 'williams',
                'fraction_internal_heating_to_surface': 1.0,
                'emissivity'                          : 0.9,
                'albedo'                              : 0.3,
                'use_real_moi'                        : True,
                'tides_on'                            : True,
                'surface_pressure'                    : 0.,
                'slices'                              : None,
                'bm_interpolation_method'             : 'mid',  # Options: mid, avg, median
                'bm_interpolation_n'                  : 100
            }
        }
    }
}
