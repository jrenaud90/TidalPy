from numba import njit

from TidalPy.exceptions import ParameterMissingError
from TidalPy.utilities.model import ModelHolder, LayerModel, ModelSearcher
from ..types import FloatArray
import numpy as np
from typing import Tuple
from . import cooling_models
from .defaults import cooling_param_defaults
import copy

# Cooling Model Finder
find_partial_melter = ModelSearcher(cooling_models, cooling_param_defaults)

class Cooling(LayerModel):

    default_config = copy.deepcopy(cooling_param_defaults)
    config_key = 'cooling'

    def __init__(self, layer):

        # Assume that gravity and density do not change after the planet is initialized.
        config = copy.deepcopy(self.default_config)
        config[layer.type]['gravity'] = layer.gravity
        config[layer.type]['density'] = layer.density
        self.default_config = config

        super().__init__(layer=layer, function_searcher=find_partial_melter, call_reinit=True)

    def _calculate(self):
        """ Calculates the cooling flux in a layer

        :return: <Tuple[ndarray's]> Heat Flux, Boundary Thickness, Rayleigh, Nusselt
        """

        temperature = self.layer.temperature
        if temperature is None:
            raise ParameterMissingError

        delta_temp = temperature - self.layer.temperature_surf


        return self.func(delta_temp, self.layer.thickness, self.layer.thermal_conductivity,
                         *self.live_inputs, *self.inputs)

    def _calculate_debug(self):

        temperature = self.layer.temperature
        if temperature is None:
            raise ParameterMissingError

        delta_temp = temperature - self.layer.temperature_surf

        assert np.all(delta_temp >= 0.)

        return self.func(delta_temp, self.layer.thickness, self.layer.thermal_conductivity,
                         *self.live_inputs, *self.inputs)