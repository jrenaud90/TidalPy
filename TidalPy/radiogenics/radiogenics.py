import copy

import numpy as np

from . import radiogenic_models
from .defaults import known_isotope_data, radiogenics_param_defaults
from .. import debug_mode
from ..exceptions import ImproperAttributeHandling, ParameterMissingError, UnknownModelError
from ..initialize import log
from ..types import float_like
from ..utilities.model import LayerModel, ModelSearcher


find_radiogenics = ModelSearcher(radiogenic_models, radiogenics_param_defaults)


class Radiogenics(LayerModel):
    default_config = copy.deepcopy(radiogenics_param_defaults)
    config_key = 'radiogenics'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=None, call_reinit=True)

        # Convert isotope information into list[tuple] format
        if 'ref_time' not in self.config:
            self.config['ref_time'] = None

        if self.model == 'isotope':
            self.isos_name = list()
            self.isos_hpr = list()
            self.isos_halflife = list()
            self.isos_massfrac = list()
            self.isos_concentration = list()

            # Isotopes may be given as a dictionary of individual isotopes or as a string pointing to one of the
            #  pre-built isotope lists.
            isotopes = self.config['isotopes']
            if type(isotopes) == str:
                if isotopes.lower() not in known_isotope_data:
                    raise UnknownModelError
                iso_datas = known_isotope_data[isotopes]
            else:
                iso_datas = isotopes

            # Different Isotope data sources my have their own reference time - extract that information
            if 'ref_time' in iso_datas:
                self.config['ref_time'] = iso_datas['ref_time']
            elif 'reference_time' in iso_datas:
                self.config['ref_time'] = iso_datas['reference_time']
            elif self.config['ref_time'] is None:
                log('No reference time provided for radiogenics, using ref_time = 0.', level='debug')
                self.config['ref_time'] = 0.

            for isotope, iso_data in iso_datas.items():
                if isotope in ['ref_time', 'reference_time']:
                    continue

                # For each isotope, extract the needed info
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

        self.func, self.inputs, self.live_input_func = find_radiogenics.find_model(self.model, self.config,
                                                                                   default_key=self.layer_type)

    def _calculate(self):
        """ Calculates the radiogenic heating of layer in which the radiogenic class is installed.

        :return: <ndarray> radiogenic heating [Watts]
        """

        return self.func(self.layer.world.time, self.layer.mass, *self.inputs)

    def _calculate_debug(self):

        time = self.layer.world.time
        assert time is not None
        assert type(time) == np.ndarray
        assert type(self.layer.mass) == float
        if self.model == 'isotope':
            if abs(time[-1] / max(self.isos_halflife)) > 1.0e4:
                log.warn('Time is much larger than radiogenic half-life - Check units of both time and half lives.')
        elif self.model == 'fixed':
            if abs(time[-1] / max(self.config['average_half_life'])) > 1.0e4:
                log.warn(
                        'Time is much larger than the fixed-average radiogenic half-life - Check units of both time and half life.')

        radio_heating = self.func(time, self.layer.mass, *self.inputs)
        if np.any(radio_heating < 0.):
            log.warn(f'Negative radiogenic heating encountered at time:\n{time[radio_heating < 0.]}')
        if np.any(radio_heating > 1.e23):
            log.warn(f'Very large radiogenic heating encountered at time:\n{time[radio_heating > 1.e23]}')

        return radio_heating

    @property
    def time(self):
        return self.layer.time

    @time.setter
    def time(self, value):
        raise ImproperAttributeHandling
