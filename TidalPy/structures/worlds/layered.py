from typing import Tuple, Dict

from .tidal import TidalWorld
from ..layers import PhysicsLayer, layers_class_by_world_class, LayerType
from ... import log
from ...exceptions import ImproperPropertyHandling, ParameterMissingError, TidalPyWorldError

BAD_LAYER_SYMBOLS = (' ', '*', '-', '/', '+', '=', '@', '#', '$', '%', '\\', '^', '&', '(', ')', '~', '`')

class LayeredWorld(TidalWorld):

    """ LayeredWorld - Construct worlds that have layers within them.

    """

    world_class = 'layered'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):

        if 'layers' not in world_config:
            log.error("Layered world's configurations do not contain layer information. "
                      "Construction can not be completed.")
            raise ParameterMissingError("Layered world's configurations do not contain layer information.")

        super().__init__(world_config, name, initialize=False)

        # Basic Configurations
        LayerClass = layers_class_by_world_class[self.world_class]
        self._layers_class = LayerClass.layer_class

        # Layer storage properties
        self._layers_by_name = dict()
        self._layers = []
        self._layer_types = []
        self._num_layers = None

        # Get layer types
        layer_i = 0
        for layer_name, layer_config in self.config['layers'].items():

            # Check for issues
            if any(bad_symbol in layer_name for bad_symbol in BAD_LAYER_SYMBOLS):
                log.error(f'An illegal symbol was found in {layer_name} for {self}.')
                raise TidalPyWorldError('Illegal symbol found in layer name.')

            # Pull out configuration info
            layer_type = layer_config['type']
            self._layer_types.append(layer_type)

            # Build Layer
            layer = LayerClass(layer_name, layer_i, self, layer_config, initialize=False)

            # Store layer in the world's containers
            self._layers.append(layer)
            self._layers_by_name[layer_name] = layer

            # Also store it in the world itself as its own attribute
            setattr(self, layer_name, layer)

        # Make layer storage immutable
        self._layers = tuple(self._layers)
        self._layer_types = tuple(self._layer_types)

        if initialize:
            self.reinit(initial_init=True, setup_simple_tides=False, reinit_layers=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, setup_simple_tides: bool = False,
               reinit_layers: bool = True):

        # Don't let parent classes initialize geometry since a LayeredWorld's mass is based on its layers' masses
        super().reinit(initial_init=initial_init, reinit_geometry=False,
                       setup_simple_tides=setup_simple_tides)

        # Pull out planet configurations
        radius = self.config['radius']
        mass = self.config.get('mass', None)

        # Layer constructor may need the planets mass and radius.
        #     So set those here (they will be reset by the set_geometry method).
        self._radius = radius
        self._mass = mass

        # Update the global tidal volume fraction
        running_tidal_fraction = 0.
        running_layer_masses = 0.

        if reinit_layers:
            # Tell the top-most layer that it is the top-most layer.
            self.layers[-1]._is_top_layer = True

            # Call reinit to the layers within this planet.
            for layer in self.layers:
                layer.reinit(initial_init)

                running_layer_masses += layer.mass
                if layer.is_tidal:
                    running_tidal_fraction += layer.tidal_scale

            self.tidal_scale = running_tidal_fraction

        if self.mass is None:
            mass = running_layer_masses
        else:
            mass = self.mass

        if reinit_geometry:
            self.set_geometry(self.radius, mass)

    def find_layer(self, layer_name: str) -> LayerType:
        """ Returns a reference to a layer with the provided name

        Layers are also stored in the planet's __dict__ and can be accessed via <world>."layer_name" as well as:
            <world>.layers_by_name (dict)
            <world>.layers (list)

        Parameters
        ----------
        layer_name : str
            Name assigned to layer in the planet's original configuration

        Returns
        -------
        layer : PhysicsLayer
            Reference to the layer class
        """

        layer = self.layers_by_name[layer_name]
        return layer

    def find_layer_by_radius(self, radius: float) -> LayerType:
        """ Returns a reference to a layer that the provided radius resides in
        If the provided radius is at the interface of two layers this method will choose the lower layer.

        Parameters
        ----------
        radius : float
            Radius [m] where the desired layer is located.

        Returns
        -------
        layer : PhysicsLayer
            Reference to the layer class at that radius.
        """

        assert radius <= self.radius

        for layer in self.layers:
            # Start from the bottom to top and look for the first layer where its radius exceeds the desired radius.
            if layer.radius >= radius:
                return layer
        raise LookupError()

    # State properties
    @property
    def layers_class(self) -> str:
        return self._layers_class

    @layers_class.setter
    def layers_class(self, value):
        raise ImproperPropertyHandling

    @property
    def layers_by_name(self) -> Dict[str, LayerType]:
        return self._layers_by_name

    @layers_by_name.setter
    def layers_by_name(self, value):
        raise ImproperPropertyHandling

    @property
    def layers(self) -> Tuple[LayerType, ...]:
        return self._layers

    @layers.setter
    def layers(self, value):
        raise ImproperPropertyHandling

    @property
    def layer_types(self) -> Tuple[str, ...]:
        return self._layer_types

    @layer_types.setter
    def layer_types(self, value):
        raise ImproperPropertyHandling

    @property
    def num_layers(self) -> int:
        return self._num_layers

    @num_layers.setter
    def num_layers(self, value):
        raise ImproperPropertyHandling


    # Dunder properties
    def __iter__(self):
        """ Planet will iterate over its layers
        Returns
        -------
        iter(self.layers)
            The iterator of the layer list.
        """

        return iter(self.layers)
