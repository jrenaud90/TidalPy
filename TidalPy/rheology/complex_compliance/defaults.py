import numpy as np

complex_compliance_defaults = {
    'ice' : {
        'model'                  : 'off',
        'voigt_compliance_offset': .2,
        'voigt_viscosity_offset' : .02,
        'alpha'                  : 1. / 3.,
        'zeta'                   : 1.,
        'critical_freq'          : 2. * np.pi / (86400. * 100.),
        'critical_freq_falloff'  : 30
        },
    'rock': {
        'model'                  : 'maxwell',
        'voigt_compliance_offset': .2,
        'voigt_viscosity_offset' : .02,
        'alpha'                  : 1. / 3.,
        'zeta'                   : 1.,
        'critical_freq'          : 2. * np.pi / (86400. * 100.),
        'critical_freq_falloff'  : 30
        },
    'iron': {
        'model'                  : 'maxwell',
        'voigt_compliance_offset': .2,
        'voigt_viscosity_offset' : .02,
        'alpha'                  : 1. / 3.,
        'zeta'                   : 1.,
        'critical_freq'          : 2. * np.pi / (86400. * 100.),
        'critical_freq_falloff'  : 30
        }
    }
