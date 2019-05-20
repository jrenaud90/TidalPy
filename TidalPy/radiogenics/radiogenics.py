from ..exceptions import ParameterMissingError
from ..utilities.classes import ModelHolder
from . import heating_models
from .. import log
from ..utilities.search import ModelSearcher
from .. import debug_mode
from ..types import float_like
import numpy as np

radiogenics_param_defaults = {
    'ice': {
        'model': 'off'
        },
    'rocky': {
            'model': 'isotope',
            'ref_time': 4600.0,
            'isotopes': {
                # Reference: Hussmann & Spohn 2004 (modified)
                'U238': {
                    'iso_mass_fraction': 0.9928,
                    'hpr': 9.48e-5,
                    'half_life': 4470.0,
                    'element_concentration': 0.012
                },
                'U235': {
                    'iso_mass_fraction': 0.0071,
                    'hpr': 5.69e-4,
                    'half_life': 704.0,
                    'element_concentration': 0.012
                },
                'Th232': {
                    'iso_mass_fraction': 0.9998,
                    'hpr': 2.69e-5,
                    'half_life': 14000.0,
                    'element_concentration': 0.04
                },
                'K40': {
                    'iso_mass_fraction': 1.19e-4,
                    'hpr': 2.92e-5,
                    'half_life': 1250.0,
                    'element_concentration': 840.0
                }
            }
        },
    'iron': {
        'model': 'off'
        }
    }

find_radiogenics = ModelSearcher(heating_models, radiogenics_param_defaults)

class Radiogenics(ModelHolder):

    name = 'Radiogenics'

    def __init__(self, layer_type: str, radiogenic_config: dict = None):

        model_searcher = ModelSearcher(heating_models, radiogenics_param_defaults[layer_type])
        model_searcher.user_parameters = radiogenic_config

        super().__init__(user_config=radiogenic_config, default_config=radiogenics_param_defaults[layer_type],
                         function_searcher=model_searcher, automate=True)

        # Convert isotope information into list[tuple] format
        self.isos_name = list()
        self.isos_hpr = list()
        self.isos_halflife = list()
        self.isos_massfrac = list()
        self.isos_concentration = list()
        for isotope, iso_data in self.config['isotopes'].items():
            self.isos_name.append(isotope)
            try:
                if debug_mode:
                    for param_name in ['hpr', 'half_life', 'iso_mass_fraction', 'element_concentration']:
                        assert type(iso_data[param_name]) in float_like
                self.isos_hpr.append(iso_data['hpr'])
                self.isos_halflife.append(iso_data['half_life'])
                self.isos_massfrac.append(iso_data['element_concentration'])
                self.isos_concentration.append(iso_data['iso_mass_fraction'])
            except KeyError:
                raise ParameterMissingError(f'One or more parameters are missing for isotope {isotope}.')
        # Once isotopes are loaded they are set and can not be appended to later. This allows the use of numba.njit
        self.isos_name = tuple(self.isos_name)
        self.isos_hpr = tuple(self.isos_hpr)
        self.isos_halflife = tuple(self.isos_halflife)
        self.isos_massfrac = tuple(self.isos_massfrac)
        self.isos_concentration = tuple(self.isos_concentration)
        self.config['iso_massfracs_of_isotope'] = self.isos_massfrac
        self.config['iso_element_concentrations'] = self.isos_concentration
        self.config['iso_halflives'] = self.isos_halflife
        self.config['iso_heat_production'] = self.isos_hpr

        self.func, self.inputs = find_radiogenics.find_model(self.model, self.config)

    def _calculate(self, time: np.ndarray, mass: float):
        """ Calculates the radiogenic heating of layer in which the radiogenic class is installed.

        :param time: <ndarray> calculation time (same units as provided half-life)
        :return: <ndarray> radiogenic heating [Watts]
        """

        return self.func(time, mass, *self.inputs)

    def _calculate_debug(self, time: np.ndarray, mass: float):

        assert type(time) == np.ndarray
        assert type(mass) in float_like
        if self.model == 'isotope':
            if abs(time[-1]/max(self.isos_halflife)) > 1.0e4:
                log.warn('Time is much larger than radiogenic half-life - Check units of both time and half lives.')
        elif self.model == 'fixed':
            if abs(time[-1]/max(self.config['average_half_life'])) > 1.0e4:
                log.warn('Time is much larger than the fixed-average radiogenic half-life - Check units of both time and half life.')

        radio_heating = self.func(time, mass, *self.inputs)
        if np.any(radio_heating < 0.):
            log.warn(f'Negative radiogenic heating encountered at time:\n{time[radio_heating < 0.]}')
        if np.any(radio_heating > 1.e23):
            log.warn(f'Very large radiogenic heating encountered at time:\n{time[radio_heating > 1.e23]}')

        return radio_heating