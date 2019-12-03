# Stored in World > Layer > 'cooling'
cooling_defaults = {
    'ice' : {
        'model'            : 'convection',
        'convection_alpha' : 1.,
        'convection_beta'  : 1. / 3.,
        'critical_rayleigh': 1600.
    },
    'rock': {
        'model'            : 'convection',
        'convection_alpha' : 1.,
        'convection_beta'  : 1. / 3.,
        'critical_rayleigh': 1100.
    },
    'iron': {
        'model': 'off'
    }
}