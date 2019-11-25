from typing import TYPE_CHECKING

import numpy as np

from ...exceptions import BadValueError
from ...performance import njit
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

        # TODO: No debug calculation has been implemented
        self._calculate_debug = self._calculate

    def _calculate(self, temperature: np.ndarray, pressure: ArrayNone = None):
        """ Wrapper for the viscosity calculator

        Parameters
        ----------
        temperature : np.ndarray
            Temperature of the layer (generally taken as the average or middle temperature) [K]
        pressure : ArrayNone = None
            Pressure of the layer (taken at the same location as temperature) [Pa]
            If set to None then a zero value for pressure will be used instead (pressure-independent viscosity)

        Returns
        -------
        viscosity : np.ndarray
            Viscosity of the layer [Pa s]
        """

        # Check if pressure was provided
        if pressure is None:
            pressure = np.asarray(0., dtype=temperature.dtype)

        viscosity = self.func(temperature, pressure, *self.inputs, *self.live_inputs)

        return viscosity


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
