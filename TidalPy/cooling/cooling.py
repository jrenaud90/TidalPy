from typing import TYPE_CHECKING

import numpy as np

from TidalPy.logger import get_logger
from TidalPy.exceptions import (IncorrectMethodToSetStateProperty, MissingAttributeError, OuterscopePropertySetError)
from TidalPy.utilities.classes.model import LayerModelHolder

from . import known_model_const_args, known_model_live_args, known_models

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray
    from TidalPy.structures.layers import PhysicalLayerType

    from .cooling_models import CoolingOutputType


log = get_logger("TidalPy")


class CoolingModel(LayerModelHolder):
    """ CoolingModel

    Cooling model provides the functionality to calculate a layer's cooling efficiency based on user provided
        parameters related to convection and conduction.

    See Also
    --------
    TidalPy.utilities.methods.model.LayerModelHolder
    """

    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = 'cooling'

    def __init__(self,
         layer: 'PhysicalLayerType', model_name: str = None, store_config_in_layer: bool = True,
        initialize: bool = True
        ):
        """ Constructor for CoolingModel class

        Parameters
        ----------
        layer : PhysicalLayerType
            The layer instance which the model should perform calculations on.
        model_name : str = None
            The user-provided model name.
        store_config_in_layer: bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `layer.config` dictionary.
        initialize : bool = True
            Determines if initial reinit should be performed on the model (loading in data from its `self.config`).
        """

        super().__init__(layer, model_name, store_config_in_layer, initialize=False)

        log.debug(f'Loading cooling model ({self.model}) into {self.layer}.')

        # State properties
        self._cooling = None
        self._cooling_flux = None
        self._boundary_layer_thickness = None
        self._rayleigh = None
        self._nusselt = None

        if initialize:
            self.reinit(initial_init=True)

    def clear_state(self):

        super().clear_state()

        self._cooling = None
        self._cooling_flux = None
        self._boundary_layer_thickness = None
        self._rayleigh = None
        self._nusselt = None

    def _calculate(self) -> 'CoolingOutputType':
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

        if self.temperature is None:
            raise MissingAttributeError(
                f"Can not calculate cooling until layer's temperature ({self.layer.name}) has been set"
                )

        if self.temperature_surf is None:
            raise MissingAttributeError(
                f"Can not calculate cooling until layer's surface temperature ({self.layer.name}) has been set. "
                f"Set the layer above its dynamic temperature, or (if top layer) set the world's surface temperature."
                )

        delta_temp = self.temperature - self.temperature_surf

        # Check if we need to use the float or array version of the complex compliance calculator
        # OPT: Put this check in the setter for the live args / freqs?
        use_float = True
        for input_ in self.live_inputs:
            if type(input_) == np.ndarray:
                use_float = False
                break
        if use_float:
            if type(delta_temp) == np.ndarray:
                use_float = False
        if use_float:
            cooling_func = self.func
        else:
            cooling_func = self.func_array

        cooling_flux, boundary_layer_thickness, rayleigh, nusselt = \
            cooling_func(delta_temp, *self.live_inputs, *self.inputs)

        self._cooling_flux = cooling_flux
        self._boundary_layer_thickness = boundary_layer_thickness
        self._rayleigh = rayleigh
        self._nusselt = nusselt
        self._cooling = self._cooling_flux * self.surface_area

        return self.cooling_flux, self.boundary_layer_thickness, self.rayleigh, self.nusselt

    def _calculate_debug(self) -> 'CoolingOutputType':

        if self.layer.surface_temperature is None:
            raise MissingAttributeError(f"Layer {self.layer.name}'s surface temperature has not been set yet.")
        if self.layer.temperature is None:
            raise MissingAttributeError(f"Layer {self.layer.name}'s average/central temperature has not been set yet.")

        return self._calculate()

    # # State properties
    @property
    def cooling(self) -> 'FloatArray':
        """ This layer's cooling [W] """
        return self._cooling

    @cooling.setter
    def cooling(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def cooling_flux(self) -> 'FloatArray':
        """ This layer's cooling flux [W m-2] """
        return self._cooling_flux

    @cooling_flux.setter
    def cooling_flux(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def boundary_layer_thickness(self) -> 'FloatArray':
        """ This layer's thermal boundary layer thickness [m] """
        return self._boundary_layer_thickness

    @boundary_layer_thickness.setter
    def boundary_layer_thickness(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def rayleigh(self) -> 'FloatArray':
        """ This layer's rayleigh number """
        return self._rayleigh

    @rayleigh.setter
    def rayleigh(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def nusselt(self) -> 'FloatArray':
        """ This layer's nusselt number """
        return self._nusselt

    @nusselt.setter
    def nusselt(self, value):
        raise IncorrectMethodToSetStateProperty

    # # Outer-scope properties
    @property
    def temperature(self):
        """ Outer-scope wrapper for layer.temperature """
        return self.layer.temperature

    @temperature.setter
    def temperature(self, value):
        raise OuterscopePropertySetError

    @property
    def temperature_surf(self):
        """ Outer-scope wrapper for layer.surface_temperature """
        return self.layer.surface_temperature

    @temperature_surf.setter
    def temperature_surf(self, value):
        raise OuterscopePropertySetError

    @property
    def surface_area(self):
        """ Outer-scope wrapper for layer.surface_area_outer """
        return self.layer.surface_area_outer

    @surface_area.setter
    def surface_area(self, value):
        raise OuterscopePropertySetError

    @property
    def viscosity(self):
        """ Outer-scope wrapper for layer.viscosity """
        return self.layer.viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise OuterscopePropertySetError

    @property
    def thermal_conductivity(self):
        """ Outer-scope wrapper for layer.thermal_conductivity """
        return self.layer.thermal_conductivity

    @thermal_conductivity.setter
    def thermal_conductivity(self, value):
        raise OuterscopePropertySetError

    @property
    def thermal_diffusivity(self):
        """ Outer-scope wrapper for layer.thermal_diffusivity """
        return self.layer.thermal_diffusivity

    @thermal_diffusivity.setter
    def thermal_diffusivity(self, value):
        raise OuterscopePropertySetError

    @property
    def thermal_expansion(self):
        """ Outer-scope wrapper for layer.thermal_expansion """
        return self.layer.thermal_expansion

    @thermal_expansion.setter
    def thermal_expansion(self, value):
        raise OuterscopePropertySetError

    @property
    def thickness(self):
        """ Outer-scope wrapper for layer.thickness """
        return self.layer.thickness

    @thickness.setter
    def thickness(self, value):
        raise OuterscopePropertySetError

    @property
    def gravity(self):
        """ Outer-scope wrapper for layer.gravity """
        return self.layer.gravity

    @gravity.setter
    def gravity(self, value):
        raise OuterscopePropertySetError

    @property
    def density_bulk(self):
        """ Outer-scope wrapper for layer.density_bulk """
        return self.layer.density_bulk

    @density_bulk.setter
    def density_bulk(self, value):
        raise OuterscopePropertySetError

    # # Aliased properties
    @property
    def blt(self):
        """ Alias for self.boundary_layer_thickness [m] """
        return self.boundary_layer_thickness

    @blt.setter
    def blt(self, value):
        self.boundary_layer_thickness = value
