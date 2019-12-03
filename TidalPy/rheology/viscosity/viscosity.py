from typing import TYPE_CHECKING

import numpy as np

from ...exceptions import ImproperAttributeHandling
from ...utilities.model import LayerModelHolder
from . import known_model_live_args, known_model_const_args, known_models
from .defaults import liquid_viscosity_defaults, solid_viscosity_defaults
from ...types import ArrayNone

if TYPE_CHECKING:
    from ...structures.layers import ThermalLayer
    from ..rheology import Rheology


class ViscosityClass(LayerModelHolder):
    """ Common Viscosity class for liquid and solid viscosity model holders.
    """

    is_liquid = False

    def __init__(self, layer: ThermalLayer, rheology_class: Rheology, model_name: str = None,
                 store_config_in_layer: bool = True):

        super().__init__(layer, model_name, store_config_in_layer)

        self.rheology_class = rheology_class

        # State properties
        self._viscosity = None

    def _calculate(self):
        """ Wrapper for the viscosity calculator

        Returns
        -------
        viscosity : np.ndarray
            Viscosity of the layer [Pa s]
        """

        viscosity = self.func(*self.live_inputs, *self.inputs)
        self._viscosity = viscosity

        return viscosity


    # State properties
    @property
    def viscosity(self):
        return self._viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise ImproperAttributeHandling

    # Outerscope references
    @property
    def temperature(self) -> np.ndarray:
        return self.layer.temperature

    @temperature.setter
    def temperature(self, value):
        raise ImproperAttributeHandling

    @property
    def pressure(self) -> np.ndarray:

        # Pressure is often turned off in various model runs, so it may not be set. If that is the case let's make sure
        #    that zeros are provided.
        _pressure = self.layer.pressure
        if _pressure is None:
            _pressure = np.asarray(0., dtype=self.temperature.dtype)

        return _pressure

    @pressure.setter
    def pressure(self, value):
        raise ImproperAttributeHandling


class LiquidViscosity(ViscosityClass):

    """ Liquid Viscosity model holder

    """

    default_config = liquid_viscosity_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('rheology', 'liquid_viscosity')
    is_liquid = True


class SolidViscosity(ViscosityClass):

    """ Solid Viscosity model holder

    """

    default_config = solid_viscosity_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('rheology', 'solid_viscosity')
    is_liquid = False
