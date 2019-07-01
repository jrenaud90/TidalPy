layer_defaults = {
    'ice' : {
        # General Layer Information
        'type'                      : 'ice',
        'slices'                    : 50,

        # Switches
        'use_surf_gravity'          : True,
        'use_tvf'                   : True,
        'is_tidal'                  : True,

        # Burnman Information
        'interp_temperature_range'  : [60., 273.15],

        # Material Information
        'material'                  : None,
        'material_source'           : None,

        # TODO: Eventually this material information should be offloaded to a material class via BM or BM Wrapper
        'stefan'                    : 0.,
        'boundary_temperature_ratio': 1.,
        'heat_fusion'               : 2.84e5,
        'shear_modulus'             : 9.2e9,
        'thermal_conductivity'      : 2.3,
        'thermal_diffusivity'       : 2.3 / (1000. * 2000.),
        'thermal_expansion'         : 5.0e-5,

    },
    'rock': {
        # General Layer Information
        'type'                      : 'rock',
        'slices'                    : 50,

        # Switches
        'use_surf_gravity'          : True,
        'use_tvf'                   : True,
        'is_tidal'                  : True,

        # Burnman Information
        'interp_temperature_range'  : [500., 2200.],

        # Material Information
        'material'                  : None,
        'material_source'           : None,
        'heat_fusion'               : None,
        'stefan'                    : 0.,
        'boundary_temperature_ratio': 1.,
        'shear_modulus'             : 6.0e10,
        'thermal_conductivity'      : 3.75,
        'thermal_diffusivity'       : 3.75 / (3250. * 1260.),
        'thermal_expansion'         : 5.2e-5,
    },
    'iron': {
        # General Layer Information
        'type'                      : 'iron',
        'slices'                    : 1,

        # Switches
        'use_surf_gravity'          : True,
        'use_tvf'                   : True,
        'is_tidal'                  : False,

        # Burnman Information
        'interp_temperature_range'  : [800., 2500.],

        # Material Information
        'material'                  : None,
        'material_source'           : None,
        'heat_fusion'               : None,
        'stefan'                    : 0.,
        'boundary_temperature_ratio': 1.1,
        'shear_modulus'             : 5.25e10,
        'thermal_conductivity'      : 7.95,
        'thermal_diffusivity'       : 7.95 / (9000. * 444.),
        'thermal_expansion'         : 1.2e-5,
    }
}

world_defaults = {
    'base' : {
        'name'                        : 'unknown_world_basetype',
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',  # Options are no_eccentricity, williams, or mendez
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'use_real_moi'                : True
    },
    'basic': {
        'name'                        : 'unknown_world_basictype',
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'use_real_moi'                : True
    },
    'star' : {
        'name'                        : 'unknown_world_startype',
        # Most of the stuff for the star class is not actually used, but having this here allows it to share the same
        # base class as the other planet types. Perhaps one day it would be a good idea to split the base classes of
        # stars and planets so that this stuff can be left off. For now though it really doesn't hurt anything.
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'use_real_moi'                : True
    },
    'tidal': {
        'name'                        : 'unknown_world_tidaltype',
        'force_spin_sync'             : True,
        'equilibrium_insolation_model': 'williams',
        'emissivity'                  : 0.9,
        'albedo'                      : 0.3,
        'quality_factor'              : 10.,
        'use_real_moi'                : True
    }
}
