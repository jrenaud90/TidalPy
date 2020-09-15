solid_viscosity_defaults = {
    'ice' : {
        'model'                     : 'arrhenius',
        # Matching Moore2006 for Volume Diffusion
        'arrhenius_coeff'           : 9.06e-8**(-1),
        'additional_temp_dependence': True,
        'stress'                    : 1.0,
        'stress_expo'               : 1.0,
        'grain_size'                : 5.0e-4,
        'grain_size_expo'           : 2.0,
        'molar_activation_energy'   : 59.4e3,
        # FIXME: I get crazy low viscosity values when I have this activation volume. For now not assuming pressure dependence.
        # 'molar_activation_volume'   : -1.3e-5
        'molar_activation_volume'   : 0.
    },
    'rock': {
        'model'                  : 'reference',
        'reference_viscosity'    : 1.0e22,
        'reference_temperature'  : 1000.0,
        'molar_activation_energy': 300000.0,
        'molar_activation_volume': 0.
    },
    'iron': {
        'model'              : 'constant',
        'reference_viscosity': 1.0e20,
    },
}

liquid_viscosity_defaults = {
    'ice' : {
        'model'                  : 'reference',
        'reference_viscosity'    : 0.89e-3,
        'reference_temperature'  : 25.0 + 273.15,
        'molar_activation_energy': 1.62e4,
        'molar_activation_volume': 0.0
    },
    'rock': {
        'model'                  : 'reference',
        'reference_viscosity'    : 0.2,
        'reference_temperature'  : 2000.0,
        'molar_activation_energy': 6.64e-20,
        'molar_activation_volume': 0.0
    },
    'iron': {
        # These values match Wijs et al 1998 (their work actually does not show much change in the liquid visc
        #    at Earth's core pressure, so a constant model may not be too incorrect).
        'model'              : 'constant',
        'reference_viscosity': 1.3e-2
    }
}
