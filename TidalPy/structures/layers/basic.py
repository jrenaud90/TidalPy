from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

from TidalPy.utilities.types import NoneType
from .defaults import layer_defaults
from ..physical import PhysicalObjSpherical
from ...exceptions import ImproperPropertyHandling, OuterscopePropertySetError, MissingArgumentError

if TYPE_CHECKING:
    from ..worlds import LayeredWorld
    from . import LayerTypes


class LayerBase(PhysicalObjSpherical):

    """ Basic Layer Class used to build more complex layer types.
    """

    default_config = layer_defaults
    layer_class = 'base'

    def __init__(self, layer_name: str, layer_index: int, world: 'LayeredWorld', layer_config: dict,
                 initialize: bool = True):

        # Load layer defaults based on layer type
        self.type = layer_config['type']
        self.default_config = self.default_config[self.type]

        # Setup Physical Layer Geometry
        super().__init__(layer_config)

        # State properties
        self._layer_index = layer_index
        self._world = world
        self._gravity = None
        self._density = None
        self._pressure = None
        self._temperature = None

        # Set up a pressure that will persist if the layer's state is cleared
        #    It is None in the case of the base layer class
        self._persistent_pressure = None

        # Other attributes
        self.name = layer_name
        self.material_name = self.config['material']
        self.tidal_scale = 1.
        self.heat_sources = None

        # Flags
        self.is_top_layer = False
        self.is_tidal = False

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):

        super().reinit()

        if not initial_init:
            self.clear_state()

        # Load in configurations
        if self.config['use_tvf']:
            self.tidal_scale = self.volume / self.world.volume
        self.is_tidal = self.config['is_tidally_active']

    def set_geometry(self, radius: float, mass: float, thickness: float = None):

        if thickness is None:
            if self.layer_index == 0:
                # Bottom layer: thickness = radius
                thickness = radius
            else:
                raise MissingArgumentError

        super().set_geometry(radius, mass, thickness)

    def clear_state(self, clear_pressure: bool = False):

        super().clear_state()

        self._temperature = None
        self._pressure = None

        if not clear_pressure:
            # Keep the old pressure so it does not have to be reinterpolated (only applicable in children classes)
            self._pressure = self._persistent_pressure

    def set_temperature(self, value):
        # The basic class implements a very simple setter. This function is mostly here to be overridden by children
        #    classes.
        self._temperature = value

    def set_pressure(self, value):
        # The basic class implements a very simple setter. This function is mostly here to be overridden by children
        #    classes.
        self._pressure = value

    # State properties
    @property
    def layer_index(self) -> int:
        return self._layer_index

    @layer_index.setter
    def layer_index(self, value):
        raise ImproperPropertyHandling('Layer index can not be changed after planet has been initialized.')

    @property
    def world(self) -> 'GeometricWorld':
        return self._world

    @world.setter
    def world(self, value):
        raise ImproperPropertyHandling('Can not change world association after a layer has been initialized.')

    @property
    def temperature(self) -> np.ndarray:
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self.set_temperature(value)

    @property
    def pressure(self) -> np.ndarray:
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        self.set_pressure(value)


    # Outer-scope properties
    # # World Class
    @property
    def layer_below(self) -> Union['LayerTypes', NoneType]:
        if self.layer_index == 0:
            return None
        else:
            return self.world.layers_by_index[self.layer_index - 1]

    @layer_below.setter
    def layer_below(self, value):
        raise OuterscopePropertySetError

    @property
    def layer_above(self) -> Union['LayerTypes', NoneType]:
        if self.layer_index == self.world.num_layers - 1:
            return None
        else:
            return self.world.layers_by_index[self.layer_index + 1]

    @layer_above.setter
    def layer_above(self, value):
        raise OuterscopePropertySetError

    @property
    def time(self):
        return self.world.time

    @time.setter
    def time(self, value):
        raise OuterscopePropertySetError


    # Dunder methods
    def __str__(self):

        if self.world is None:
            text = f'[Layer {self.name.title()}:{self.type.title()} no world]'
        else:
            text = f'[Layer {self.name.title()}:{self.type.title()} in {self.world}]'
        return text

    def __repr__(self):

        text = f'{self.__class__} object at {hex(id(self))}'
        if 'name' in self.__dict__:
            if self.name is not None:
                text = f'{self.name} ' + text
        if self.world is not None:
            text += f'; stored in {repr(self.world)}'
        else:
            text += '; not associated with a world'

        return text
