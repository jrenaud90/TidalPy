layer_defaults = {
    'ice' : {
        # General Layer Information
        'type'                         : 'ice',

        # Switches
        'use_surface_gravity'          : True,
        'use_bulk_density'             : True,
        'use_tidal_vol_frac'           : True,
        'is_tidally_active'            : True,
        'use_pressure_in_strength_calc': False,

        # Burnman Information
        'interp_temperature_range'     : [20., 273.15],

        # Material Information
        'material'                     : None,
        'material_source'              : None,
        'slices'                       : 40,

        # TODO: Eventually this material information should be offloaded to a material class via BM or BM Wrapper
        'stefan'                       : 0.,
        'boundary_temperature_ratio'   : 1.,
        'heat_fusion'                  : 2.84e5,
        'shear_modulus'                : 9.2e9,
        'thermal_conductivity'         : 2.3,
        'thermal_diffusivity'          : 2.3 / (1000. * 2000.),
        'thermal_expansion'            : 5.0e-5,

    },
    'rock': {
        # General Layer Information
        'type'                         : 'rock',

        # Switches
        'use_surface_gravity'          : True,
        'use_bulk_density'             : True,
        'use_tidal_vol_frac'           : True,
        'is_tidally_active'            : True,
        'use_pressure_in_strength_calc': False,
        'slices'                       : 40,

        # Burnman Information
        'interp_temperature_range'     : [100., 3200.],

        # Material Information
        'material'                     : None,
        'material_source'              : None,
        'heat_fusion'                  : None,
        'stefan'                       : 0.,
        'boundary_temperature_ratio'   : 1.,
        'shear_modulus'                : 6.0e10,
        'thermal_conductivity'         : 3.75,
        'thermal_diffusivity'          : 3.75 / (3250. * 1260.),
        'thermal_expansion'            : 5.2e-5,
    },
    'iron': {
        # General Layer Information
        'type'                         : 'iron',

        # Switches
        'use_surface_gravity'          : True,
        'use_bulk_density'             : True,
        'use_tidal_vol_frac'           : True,
        'is_tidally_active'            : False,
        'use_pressure_in_strength_calc': False,
        'slices'                       : 40,

        # Burnman Information
        'interp_temperature_range'     : [500., 5000.],

        # Material Information
        'material'                     : None,
        'material_source'              : None,
        'heat_fusion'                  : None,
        'stefan'                       : 0.,
        'boundary_temperature_ratio'   : 1.1,
        'shear_modulus'                : 5.25e10,
        'thermal_conductivity'         : 7.95,
        'thermal_diffusivity'          : 7.95 / (9000. * 444.),
        'thermal_expansion'            : 1.2e-5,
    }
}
