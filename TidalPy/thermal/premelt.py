from ..utilities.search import ModelSearcher
from . import viscosity

viscosity_param_defaults = {
    'ice': {
        'boundary_diffusion':
            {'arrhenius_coeff': }
arrhenius_coeff, stress, stress_expo, grain_size, grain_size_expo, molar_activation_energy, molar_activation_volume
        },

    'rocky': {

        },

    'iron': {

        }
    }

find_viscosity = ModelSearcher(viscosity, viscosity_param_defaults)