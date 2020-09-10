from .tidal import TidalWorld
from ...exceptions import ImproperPropertyHandling
from ..layers import PhysicsLayer

class LayeredWorld(TidalWorld):

    """ LayeredWorld - Construct worlds that have layers within them.

    """

    world_class = 'layered'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):

        super().__init__(self, world_config, name, initialize=False)

        # Layer storage properties
        self._layers_by_name = dict()
        self._layers = None

    def find_layer(self, layer_name: str) -> PhysicsLayer:
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

    def find_layer_by_radius(self, radius: float) -> PhysicsLayer:
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
    def layers_by_name(self):
        return self._layers_by_name

    @layers_by_name.setter
    def layers_by_name(self, value):
        raise ImproperPropertyHandling

    @property
    def layers_by_order(self):
        return self._layers_by_order

    @layers_by_order.setter
    def layers_by_order(self, value):
        raise ImproperPropertyHandling

    @property
    def layers(self):
        return self._layers

    @layers.setter
    def layers(self, value):
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




x={"layers": {
        "Core": {
            "type": "iron",
            "is_tidal": false,
            "radius": 810.0e3,
            "material": [
                "Pyrite",
                "Fe_Dewaele"
            ],
            "material_source": [
                "TidalPy",
                "other"
            ],
            "material_fractions": [
                0.5,
                0.5
            ],
            "temperature_mode": "user-defined",
            "temperature_fixed": 1800.0
        },
        "Mantle": {
            "type": "rock",
            "is_tidal": true,
            "radius": 1821.49e3,
            "material": [
                "forsterite",
                "mg_perovskite"
            ],
            "material_source": [
                "SLB_2011",
                "SLB_2011"
            ],
            "material_fractions": [
                0.65,
                0.35
            ],
            "temperature_mode": "adiabatic",
            "temperature_top": 1800.0,
            "surface_temperature": 100.0
        }}