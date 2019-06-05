from TidalPy import debug_mode
from TidalPy.exceptions import ParameterMissingError, ReinitError, ImproperAttributeHandling
from TidalPy.structures.layers import construct_layer
from .physical import PhysicalObjSpherical
from ..physics.stellar import luminosity_from_mass, efftemp_from_luminosity, luminosity_from_efftemp
from ..graphics.planet_plot import geotherm_plot
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

    def __init__(self, planet_config: dict, burnman_world: burnman.Planet, bm_layers: list, automate: bool = True):

        super().__init__(planet_config, automate=automate)

        self.bm_world = burnman_world
        self.bm_layers = bm_layers

        # Setup geometry based on burnman
        self.set_geometry(self.bm_world.radius_planet, self.bm_world.mass)

        # Setup layers
        self.layers_byname = dict()
        self.layers = list()
        for bm_layer in self.bm_layers:
            layer_name = bm_layer.name
            layer_config = self.config['layers'][layer_name]
            layer = construct_layer(layer_name, bm_layer, layer_config)
            self.layers_byname[layer_name] = layer
            self.layers.append(layer)
            setattr(self, layer.name.lower(), layer)

        # Provide references to other layers
        for layer_below, layer, layer_above in zip([None] + self.layers[:-1], self.layers, self.layers[1:] + [None]):
            layer.layer_below = layer_below
            layer.layer_above = layer_above

        # Concatenate the slice information from the layers
        gravity_list = [layer.gravity_slices for layer in self]
        pressure_list = [layer.pressure_slices for layer in self]
        density_list = [layer.density_slices for layer in self]
        radius_list = [layer.radii for layer in self]
        self.gravity_slices = np.concatenate(gravity_list)
        self.pressure_slices = np.concatenate(pressure_list)
        self.density_slices = np.concatenate(density_list)
        self.radii = np.concatenate(radius_list)

        # Independent State variables
        self._spin_freq = None
        self._orbital_freq = None
        self._time = None
        self._eccentricity = None
        self._inclination = None

        # Dependent variables
        self._tidal_modes = None
        self._tidal_coeff_heating = None
        self._tidal_coeff_torque = None


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

    def get_layer(self, layer_name: str):
        """ Returns a layer with a given name

        :param layer_name: <str> layer's name"""

        return self.layers_byname[layer_name]

    def get_layer_by_radius(self, radius: float):
        """ Returns a layer in which this radius lies

        :param radius: <float> radius inside the world
        """

        assert radius <= self.radius

        for layer in self.layers:
            if layer.radius >= radius:
                return layer
        raise LookupError()

    def update_modes(self):

        self._tidal_modes = tuple(self.orbital_freq)
        self._tidal_coeff_heating = tuple(7. * self.eccentricity )
        self._tidal_coeff_torque = tuple()

        for layer in self:
            layer.tidal_modes = self._tidal_modes

    def paint(self, depth_plot: bool = False, auto_show: bool = False):
        """ Create a geotherm or depth plot of the planet's gravity, pressure, and density

        :param depth_plot: <bool> (=False) If true then plot will be vs. depth instead of vs. radius
        :param auto_show:  <bool> (=False) If true then plt.show will be called
        :return: matplotlib figure
        """

        return geotherm_plot(self.radii, self.gravity_slices, self.pressure_slices, self.density_slices,
                             bulk_density=self.density_bulk,
                             planet_radius=self.radius, depth_plot=depth_plot, auto_show=auto_show)

    @property
    def time(self):
        return self.time

    @time.setter
    def time(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._time = value
        for layer in self:
            layer.time = value

    @property
    def spin_freq(self):
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._spin_freq = value

    @property
    def orbital_freq(self):
        return self._orbital_freq

    @orbital_freq.setter
    def orbital_freq(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._orbital_freq = value

    @property
    def eccentricity(self):
        return self._eccentricity

    @eccentricity.setter
    def eccentricity(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._eccentricity = value

    @property
    def inclination(self):
        return self._orbital_freq

    @inclination.setter
    def inclination(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._inclination = value

    @property
    def tidal_modes(self):
        return self._tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_coeff_heating(self):
        return self._tidal_coeff_heating

    @tidal_coeff_heating.setter
    def tidal_coeff_heating(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_coeff_torque(self):
        return self._tidal_coeff_torque

    @tidal_coeff_torque.setter
    def tidal_coeff_torque(self, value):
        raise ImproperAttributeHandling

    def __iter__(self):

        return iter(self.layers)


class MultiModeWorld(BurnManWorld):

    def update_modes(self):

        self._tidal_modes = tuple(self.orbital_freq)
        self._tidal_coeff_heating = tuple()
        self._tidal_coeff_torque = tuple()

        for layer in self:
            layer.tidal_modes = self._tidal_modes


world_types = {
    'star': Star,
    'basic': BasicWorld,
    'simple': BasicWorld,
    'gas': BasicWorld,
    'burnman': BurnManWorld,
    'tidal': BurnManWorld
}