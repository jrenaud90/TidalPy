from typing import TYPE_CHECKING, Tuple

import numpy as np

from ..exceptions import MissingAttributeError
from ..utilities.model import LayerModelHolder
from . import known_models, known_model_live_args, known_model_const_args
from .defaults import cooling_defaults

if TYPE_CHECKING:
    from ..structures.layers import ThermalLayer


class CoolingModel(LayerModelHolder):

    """ Cooling Model Class - Child of LayerModelHolder Class

    Cooling model provides the functionality to calculate a layer's cooling efficiency based on user provided
        parameters related to convection and conduction.
    """

    default_config = cooling_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('cooling', )

    def __init__(self, layer: ThermalLayer, model_name: str = None, store_config_in_layer: bool = True):

        super().__init__(layer, model_name, store_config_in_layer)

    def _calculate(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """ Calculate layer cooling based on the layer's state properties.

        Returns
        -------
        cooling_flux : np.ndarray
            Heat flux leaving the layer [W m-2]
        boundary_layer_thickness : np.ndarray
            Thickness of boundary layer (if any) [m]
        rayleigh : np.ndarray
            Rayleigh number (only non-zero for convection)
        nusselt : np.ndarray
            Nusselt number (only != 1 for convection)
        """

        surface_temp = self.layer.temperature_surf
        temp = self.layer.temperature

        delta_temp = temp - surface_temp

        return self.func(delta_temp)

    def _calculate_debug(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

        if self.layer.temperature_surf is None:
            raise MissingAttributeError(f"Layer {self.layer.name}'s surface temperature has not been set yet.")
        if self.layer.temperature is None:
            raise MissingAttributeError(f"Layer {self.layer.name}'s average/central temperature has not been set yet.")

        return self._calculate()