radiogenics_param_defaults = {
    'ice' : {
        'model': 'off'
    },
    'rock': {
        'model'   : 'isotope',
        'ref_time': 4600.0,
        'isotopes': {
            # Reference: Hussmann & Spohn 2004 (modified)
            'U238' : {
                'iso_mass_fraction'    : 0.9928,
                'hpr'                  : 9.48e-5,
                'half_life'            : 4470.0,
                'element_concentration': 0.012
            },
            'U235' : {
                'iso_mass_fraction'    : 0.0071,
                'hpr'                  : 5.69e-4,
                'half_life'            : 704.0,
                'element_concentration': 0.012
            },
            'Th232': {
                'iso_mass_fraction'    : 0.9998,
                'hpr'                  : 2.69e-5,
                'half_life'            : 14000.0,
                'element_concentration': 0.04
            },
            'K40'  : {
                'iso_mass_fraction'    : 1.19e-4,
                'hpr'                  : 2.92e-5,
                'half_life'            : 1250.0,
                'element_concentration': 840.0
            }
        }
    },
    'iron': {
        'model': 'off'
    }
}

known_iso_data = radiogenics_param_defaults['rock']['isotopes']

standard_isotope_input = (
    tuple([iso_data['iso_mass_fraction'] for _, iso_data in known_iso_data.items()]),
    tuple([iso_data['element_concentration'] for _, iso_data in known_iso_data.items()]),
    tuple([iso_data['half_life'] for _, iso_data in known_iso_data.items()]),
    tuple([iso_data['hpr'] for _, iso_data in known_iso_data.items()]),
    4600.
)
