from TidalPy import debug_mode
from TidalPy.exceptions import ParameterMissingError, ReinitError
from TidalPy.structures.layers import construct_layer
from .physical import PhysicalObjSpherical
from ..physics.stellar import luminosity_from_mass, efftemp_from_luminosity, luminosity_from_efftemp

import burnman
import numpy as np

class WorldBase(PhysicalObjSpherical):

    def __init__(self, world_config: dict, automate: bool = True):

        super().__init__(config=world_config, automate=automate)


    def set_geometry(self, radius: float, mass: float):

        super().set_geometry(radius, mass, thickness=radius)


class BasicWorld(WorldBase):

    def __init__(self, planet_config: dict, automate: bool = True):
        super().__init__(planet_config, automate=automate)

        self.set_geometry(planet_config['radius'], planet_config['mass'])


class Star(BasicWorld):


    def __init__(self, star_config: dict, automate: bool = True):

        super().__init__(star_config, automate=automate)

        # Star-specific attributes
        self.luminosity = None
        self.effective_temperature = None

    def init(self):

        self.luminosity = self.config.get('luminosity', None)
        self.effective_temperature = self.config.get('effective_temperature', None)

        if self.luminosity is None:
            # If no luminosity provided: Try to convert effective surface temperature
            if self.effective_temperature is None:
                # if that fails, try to estimate from mass
                self.luminosity = luminosity_from_mass(self.mass)
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)
            else:
                self.luminosity = luminosity_from_efftemp(self.effective_temperature, self.radius)
        else:
            raise ParameterMissingError


class BurnManWorld(WorldBase):

    def __init__(self, planet_config: dict, burnman_world: burnman.Planet, bm_layers: dict,
                 automate: bool = False):

        super().__init__(planet_config, automate=automate)

        self.bm_world = burnman_world
        self.bm_layers = bm_layers

        # Setup geometry based on burnman
        self.set_geometry(self.bm_world.radius_planet, self.bm_world.mass)

        # Setup layers
        self.layers_byname = dict()
        self.layers = list()
        mass_below = 0.
        for layer_name, bm_layer in self.bm_layers:
            layer_config = self.config['layers'][layer_name]
            layer = construct_layer(layer_name, bm_layer, mass_below, layer_config)
            self.layers_byname[layer_name] = layer
            self.layers.append(layer)
            mass_below += layer.mass

        # Provide references to other layers
        for layer_below, layer, layer_above in zip([None] + self.layers[:-1], self.layers, self.layers[1:] + [None]):
            layer.layer_below = layer_below
            layer.layer_above = layer_above


    def reinit(self):

        # Determine if planet must be made again due to config changes that impact burnman results
        for layer_name, layer_dict in self.config['layers'].items():
            old_layer_dict = self._old_config[layer_name]
            for critical_attribute in ['material', 'material_source', 'type']:
                if layer_dict[critical_attribute] != old_layer_dict[critical_attribute]:
                    raise ReinitError

        # Load in new configs into layers and call their reinits
        for layer in self:
            layer.config = self.config['layers'][layer.name]
            layer.reinit()

        super().reinit()

    @property
    def time(self):
        return self.time

    @time.setter
    def time(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._time = value
        for layer in self:
            layer._time = value

    @property
    def spin_freq(self):
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._spin_freq = value
        for layer in self:
            layer._spin_freq = value

    @property
    def orbital_freq(self):
        return self._orbital_freq

    @orbital_freq.setter
    def orbital_freq(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._orbital_freq = value
        for layer in self:
            layer._orbital_freq = value

    def __iter__(self):

        return iter(self.layers)


world_types = {
    'star': Star,
    'basic': BasicWorld,
    'simple': BasicWorld,
    'gas': BasicWorld,
    'burnman': BurnManWorld,
    'tidal': BurnManWorld
}