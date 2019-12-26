from typing import TYPE_CHECKING

import numpy as np

from .defaults import known_isotope_data, radiogenics_defaults
from . import known_models, known_model_const_args, known_model_live_args
from .. import debug_mode
from ..exceptions import (ImproperAttributeHandling, ParameterMissingError, UnknownModelError, AttributeNotSetError,
                          IncorrectAttributeType)
from ..initialize import log
from ..types import float_like
from ..utilities.model import LayerModelHolder


class Radiogenics(LayerModelHolder):

    """ Radiogenic Model Class - Child of LayerModelHolder Class

    Radiogenic model provides the functionality to calculate a layer's heating due to radioactive isotopes based on
        user provided parameters related to convection and conduction.
    """

    default_config = radiogenics_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = 'radiogenics'

    def __init__(self, layer, model_name: str = None, store_config_in_layer: bool = True):

        # Set auto_build_inputs to False so that the functions can be built at the end of this __init__ once a few
        #    more parameters are loaded into the class's config.
        super().__init__(layer, model_name, store_config_in_layer, auto_build_inputs=False)

        # State attributes
        self._heating = None

        log(f'Loading Radiogenics into {self.layer}...', level='info')

        # Convert isotope information into list[tuple] format
        if 'ref_time' not in self.config:
            if 'reference_time' in self.config:
                self.config['ref_time'] = self.config['reference_time']
            else:
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

            # Different Isotope data sources may have their own reference time - extract that information
            new_ref_time = None
            if 'ref_time' in iso_datas:
                new_ref_time = iso_datas['ref_time']
            elif 'reference_time' in iso_datas:
                new_ref_time = iso_datas['reference_time']
            elif self.config['ref_time'] is None:
                log('\tNo reference time provided for radiogenics, using ref_time = 0.', level='debug')
                self.config['ref_time'] = 0.

            if new_ref_time is not None:
                if self.config['ref_time'] is not None:
                    log('\tRadiogenic reference time overwritten by isotope reference time.', level='debug')
                self.config['ref_time'] = new_ref_time

            # Build isotope inputs and store them in radiogenics.config - these will end up in the self.inputs used in
            #    self.calculate
            for isotope, iso_data in iso_datas.items():
                if isotope in ['ref_time', 'reference_time']:
                    # Skip non-isotope keys
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

            # Once isotopes are loaded they are set and can not be appended or changed to later.
            #    This allows the use of numba.njit.
            #    To emphasize this we will convert these lists into tuples for final storage.
            self.isos_name = tuple(self.isos_name)
            self.isos_hpr = tuple(self.isos_hpr)
            self.isos_halflife = tuple(self.isos_halflife)
            self.isos_massfrac = tuple(self.isos_massfrac)
            self.isos_concentration = tuple(self.isos_concentration)
            self.config['iso_massfracs_of_isotope'] = self.isos_massfrac
            self.config['iso_element_concentrations'] = self.isos_concentration
            self.config['iso_halflives'] = self.isos_halflife
            self.config['iso_heat_production'] = self.isos_hpr

        # Finish loading the model functions and inputs
        self.build_inputs()

    def _calculate(self) -> np.ndarray:
        """ Calculates the radiogenic heating of layer in which the radiogenic class is installed.

        Returns
        -------
        radiogenic_heating : np.ndarray
            Radiogenic heating [W]
        """

        radiogenic_heating = self.func(*self.live_inputs, *self.inputs)
        self._heating = radiogenic_heating

        return radiogenic_heating

    def _calculate_debug(self):

        # Attribute checks
        if self.time is None:
            raise AttributeNotSetError
        if self.mass is None:
            raise AttributeNotSetError
        if type(self.time) != np.ndarray:
            raise IncorrectAttributeType
        if type(self.mass) != float:
            raise IncorrectAttributeType

        # Value checks
        if self.model == 'isotope':
            if abs(self.time[-1] / max(self.isos_halflife)) > 1.0e4:
                log.warn('Time is much larger than radiogenic half-life - Check units of the time array and half lives.')
        elif self.model == 'fixed':
            if abs(self.time[-1] / max(self.config['average_half_life'])) > 1.0e4:
                log.warn('Time is much larger than the fixed-average radiogenic half-life - '
                         'Check units of the time array and fixed half-life.')

        # Calculate and perform more value checks
        radio_heating = self.func(self.time, self.layer.mass, *self.inputs)
        if np.any(radio_heating < 0.):
            log.warn(f'Negative radiogenic heating encountered at time:\n{self.time[radio_heating < 0.]}')
        if np.any(radio_heating > 1.e23):
            log.warn(f'Very large radiogenic heating encountered at time:\n{self.time[radio_heating > 1.e23]}')

        return radio_heating

    # State properties
    @property
    def heating(self) -> np.ndarray:
        return self._heating

    @heating.setter
    def heating(self, value):
        raise ImproperAttributeHandling

    # Outerscope properties
    @property
    def time(self):
        return self.layer.world.time

    @time.setter
    def time(self, value):
        raise ImproperAttributeHandling

    @property
    def mass(self):
        return self.layer.radiogenic_mass

    @mass.setter
    def mass(self, value):
        raise ImproperAttributeHandling
