import copy

from ..constants import ppm

known_isotope_data = {
    'modern_day_chondritic': {
        'ref_time': 4600.0,
        # Based off Hussmann & Spohn 2004 and Turcotte & Schubert 2001
        'U238'    : {
            'iso_mass_fraction'    : 0.9928,
            'hpr'                  : 9.48e-5,
            'half_life'            : 4470.0,
            'element_concentration': 0.012 * ppm
        },
        'U235'    : {
            'iso_mass_fraction'    : 0.0071,
            'hpr'                  : 5.69e-4,
            'half_life'            : 704.0,
            'element_concentration': 0.012 * ppm
        },
        'Th232'   : {
            'iso_mass_fraction'    : 0.9998,
            'hpr'                  : 2.69e-5,
            'half_life'            : 14000.0,
            'element_concentration': 0.04 * ppm
        },
        'K40'     : {
            'iso_mass_fraction'    : 1.19e-4,
            'hpr'                  : 2.92e-5,
            'half_life'            : 1250.0,
            'element_concentration': 840.0 * ppm
        }
    },
    'LLRI_and_SLRI'        : {
        'ref_time': 4600.0,
        # Based off Castillo-Rogez et al 2007
        'U238'    : {
            'iso_mass_fraction'    : 0.9928,
            'hpr'                  : 94.65e-6,
            'half_life'            : 4468.0,
            'element_concentration': 0.026 * ppm
        },
        'U235'    : {
            'iso_mass_fraction'    : 0.0071,
            'hpr'                  : 568.7e-6,
            'half_life'            : 703.81,
            'element_concentration': 0.0082 * ppm
        },
        'Th232'   : {
            'iso_mass_fraction'    : 1.0,
            'hpr'                  : 26.38e-6,
            'half_life'            : 14025.0,
            'element_concentration': 0.0538 * ppm
        },
        'K40'     : {
            'iso_mass_fraction'    : 1.176e-4,
            'hpr'                  : 29.17e-6,
            'half_life'            : 1277.0,
            'element_concentration': 1.104 * ppm
        },
        'Mn53'    : {
            'iso_mass_fraction'    : 2.0e-5,
            'hpr'                  : 0.027,
            'half_life'            : 3.7,
            'element_concentration': 0.0257 * ppm
        },
        'Fe60'    : {
            'iso_mass_fraction'    : 1.0e-6,
            'hpr'                  : 0.07,
            'half_life'            : 1.5,
            'element_concentration': 0.1 * ppm
        },
        'Al26'    : {
            'iso_mass_fraction'    : 5.0e-5,
            'hpr'                  : 0.146,
            'half_life'            : 0.72,
            'element_concentration': 0.6 * ppm
        }
    }
}

radiogenics_defaults = {
    'ice' : {
        'model'                              : 'off',
        'use_full_layer_mass_for_radiogenics': True,
        'radiogenic_layer_mass_fraction'     : 1.
    },
    'rock': {
        'model'                         : 'isotope',
        # isotopes can either be a string that references one of the pre-defined sets above,
        #  or a custom set (use the above format).
        'isotopes'                      : 'modern_day_chondritic',
        'radiogenic_layer_mass_fraction': 1.
    },
    'iron': {
        'model'                              : 'off',
        'use_full_layer_mass_for_radiogenics': True,
        'radiogenic_layer_mass_fraction'     : 1.
    }
}

standard_isotope_data = copy.deepcopy(known_isotope_data['modern_day_chondritic'])
del standard_isotope_data['ref_time']

standard_isotope_input = (
    tuple([iso_data['iso_mass_fraction'] for _, iso_data in standard_isotope_data.items()]),
    tuple([iso_data['element_concentration'] for _, iso_data in standard_isotope_data.items()]),
    tuple([iso_data['half_life'] for _, iso_data in standard_isotope_data.items()]),
    tuple([iso_data['hpr'] for _, iso_data in standard_isotope_data.items()]),
    4600.
)
