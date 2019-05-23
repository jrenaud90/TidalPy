from TidalPy.exceptions import ParameterMissingError
from .physical import PhysicalObjSpherical
from ..physics.stellar import luminosity_from_mass, efftemp_from_luminosity, luminosity_from_efftemp

import burnman

class WorldBase(PhysicalObjSpherical):

    def __init__(self, world_config: dict, automate: bool = True):

        super().__init__(config=world_config, automate=automate)


    def set_geometry(self, radius: float, mass: float):

        super().set_geometry(radius, mass, thickness=radius)


class BasicWorld(WorldBase):

    def __init__(self, planet_config: dict):
        super().__init__(planet_config)

        self.set_geometry(planet_config['radius'], planet_config['mass'])


class Star(BasicWorld):


    def __init__(self, star_config: dict):

        super().__init__(star_config)

        # Star-specific attributes
        self.luminosity = self.config.get('luminosity', None)
        self.effective_temperature = self.config.get('effective_temperature', None)

    def set_luminosity(self, luminosity: float = None):
        """ Sets the star's luminosity if value is provided, otherwise an estimation will be made

        :param luminosity: <float> Star's luminosity [Watts]
        :return:           <float> Star's luminosity [Watts]
        """

        eff_temp = self.config.get('effective_temperature', None)

        if luminosity is None:
            # If no luminosity provided: Try to convert effective surface temperature
            if self.effective_temperature is not None:
                luminosity = luminosity_from_efftemp(eff_temp, self.radius)
            else:
                # if that fails, try to estimate from mass
                luminosity = luminosity_from_mass(self.mass)

        # Updating luminosity will change the effective temperature
        self.luminosity = luminosity
        self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)

        return self.luminosity


class BurnManWorld(WorldBase):

    def __init__(self, planet_config: dict, burnman_world: burnman.Planet, bm_layers: dict,
                 automate: bool = False):

        super().__init__(planet_config, automate=automate)

        self.bm_world = burnman_world
        self.bm_layers = bm_layers

        # Setup geometry based on burnman
        self.set_geometry(self.bm_world.radius_planet, self.bm_world.mass)

        # Setup layers
        self.layers = dict()
        self.layers_byorder = list()
        for layer_name, bm_layer in self.bm_layers:
            self.layers[layer_name] = Layer()

    def reinit(self):
        # Determine if planet must be made again due to config changes that impact burnman results

        # Load in new configs into layers and call their reinits