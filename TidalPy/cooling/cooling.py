from typing import TYPE_CHECKING, Tuple

import numpy as np

from ..initialize import log
from ..exceptions import MissingAttributeError, ImproperAttributeHandling, OuterscopeAttributeSetError
from ..utilities.model import LayerModelHolder
from . import known_models, known_model_live_args, known_model_const_args
from .defaults import cooling_defaults


class CoolingModel(LayerModelHolder):

    """ Cooling Model Class - Child of LayerModelHolder Class

    Cooling model provides the functionality to calculate a layer's cooling efficiency based on user provided
        parameters related to convection and conduction.
    """

    default_config = cooling_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = 'cooling'

    def __init__(self, layer, model_name: str = None, store_config_in_layer: bool = True):

        super().__init__(layer, model_name, store_config_in_layer)

        # State properties
        self._cooling = None
        self._cooling_flux = None
        self._boundary_layer_thickness = None
        self._rayleigh = None
        self._nusselt = None

        # Report model building
        log(f'Cooling model build in {self.layer}: {self.model}.')

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

        cooling_flux, boundary_layer_thickness, rayleigh, nusselt = \
            self.func(delta_temp, *self.live_inputs, *self.inputs)

        self._cooling_flux = cooling_flux
        self._boundary_layer_thickness = boundary_layer_thickness
        self._rayleigh = rayleigh
        self._nusselt = nusselt
        self._cooling = self._cooling_flux * self.layer.surface_area_outer

        return self.cooling_flux, self.boundary_layer_thickness, self.rayleigh, self.nusselt

    def _calculate_debug(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

        if self.layer.temperature_surf is None:
            raise MissingAttributeError(f"Layer {self.layer.name}'s surface temperature has not been set yet.")
        if self.layer.temperature is None:
            raise MissingAttributeError(f"Layer {self.layer.name}'s average/central temperature has not been set yet.")

        return self._calculate()

    # State properties
    @property
    def cooling(self) -> np.ndarray:
        return self._cooling

    @cooling.setter
    def cooling(self, value):
        raise ImproperAttributeHandling

    @property
    def cooling_flux(self) -> np.ndarray:
        return self._cooling_flux

    @cooling_flux.setter
    def cooling_flux(self, value):
        raise ImproperAttributeHandling

    @property
    def boundary_layer_thickness(self) -> np.ndarray:
        return self._boundary_layer_thickness

    @boundary_layer_thickness.setter
    def boundary_layer_thickness(self, value):
        raise ImproperAttributeHandling

    @property
    def rayleigh(self) -> np.ndarray:
        return self._rayleigh

    @rayleigh.setter
    def rayleigh(self, value):
        raise ImproperAttributeHandling

    @property
    def nusselt(self) -> np.ndarray:
        return self._nusselt

    @nusselt.setter
    def nusselt(self, value):
        raise ImproperAttributeHandling

    # Outerscope properties
    @property
    def viscosity(self):
        return self.layer.viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def thermal_conductivity(self):
        return self.layer.thermal_conductivity

    @thermal_conductivity.setter
    def thermal_conductivity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def thermal_diffusivity(self):
        return self.layer.thermal_diffusivity

    @thermal_diffusivity.setter
    def thermal_diffusivity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def thermal_expansion(self):
        return self.layer.thermal_expansion

    @thermal_expansion.setter
    def thermal_expansion(self, value):
        raise OuterscopeAttributeSetError

    @property
    def thickness(self):
        return self.layer.thickness

    @thickness.setter
    def thickness(self, value):
        raise OuterscopeAttributeSetError

    @property
    def gravity(self):
        # TODO: Should this be surface or central gravity?
        return self.layer.gravity_surface

    @gravity.setter
    def gravity(self, value):
        raise OuterscopeAttributeSetError

    @property
    def density_bulk(self):
        return self.layer.density_bulk

    @density_bulk.setter
    def density_bulk(self, value):
        raise OuterscopeAttributeSetError

    # Alias properties
    @property
    def blt(self):
        return self.boundary_layer_thickness

    @blt.setter
    def blt(self, value):
        self.boundary_layer_thickness = value
