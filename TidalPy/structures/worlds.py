from __future__ import annotations
from typing import Tuple

from TidalPy import debug_mode
from TidalPy.exceptions import ParameterMissingError, ReinitError, ImproperAttributeHandling, BadArrayShape
from TidalPy.orbit.modes import spin_sync_modes, nsr_modes
from TidalPy.structures.layers import construct_layer, TidalLayer
from .physical import PhysicalObjSpherical
from ..physics.stellar import luminosity_from_mass, efftemp_from_luminosity, luminosity_from_efftemp
from ..graphics.planet_plot import geotherm_plot
import burnman
import numpy as np
from .defaults import world_defaults

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..orbit.orbit import OrbitBase

class WorldBase(PhysicalObjSpherical):

    type = 'base'

    def __init__(self, world_config: dict, automate: bool = True):

        # Load in defaults
        self.default_config = world_defaults[self.type]

        super().__init__(config=world_config, automate=automate)

        self.pyname = f'{self.__class__}_{self.type}'

        # Pull out switch information
        self.is_spin_sync = self.config['force_spin_sync']

        # Independent State variables
        self._spin_freq = None
        self._time = None

        # Orbit reference
        self._orbit = None  # type: OrbitBase
        self.orbit_location = None


    def set_geometry(self, radius: float, mass: float):

        super().set_geometry(radius, mass, thickness=radius)

    def orbit_change(self):
        """ Performs state updates whenever the planet's orbital parameters are changed

        """
        # Base class has nothing to update.
        pass

    def set_state(self, orbital_freq: np.ndarray = None, semi_major_axis: np.ndarray = None,
                  eccentricity: np.ndarray = None, inclination: np.ndarray = None,
                  spin_freq: np.ndarray = None):
        """ Set multiple orbital parameters at once, this reduces the number of calls to self.orbit_change

        :param orbital_freq:    <ndarray> (Optional, exclusive w/ semi_a)
        :param semi_major_axis: <ndarray> (Optional, exclusive w/ orb_freq)
        :param eccentricity:    <ndarray> (Optional)
        :param inclination:     <ndarray> (Optional)
        :param spin_freq:       <ndarray> (Optional)
        """

        self.orbit.set_state(self.orbit_location, orbital_freq, semi_major_axis, eccentricity, inclination,
                             set_by_planet=True)
        if spin_freq is not None:
            if type(spin_freq) != np.ndarray:
                spin_freq = np.asarray(spin_freq)
            self._spin_freq = spin_freq

        self.orbit_change()

    @property
    def orbit(self):
        if debug_mode and self._orbit is None:
            raise ParameterMissingError
        return self._orbit

    @orbit.setter
    def orbit(self, value):
        raise ImproperAttributeHandling('Orbit is automatically passed to a planet when a new orbit class is created.')

    # Independent state variables
    @property
    def time(self):
        return self.time

    @time.setter
    def time(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._time = value

    @property
    def spin_freq(self):
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._spin_freq = value
        if not self.is_spin_sync:
            # NSR will change tidal modes
            self.orbit_change()

    # Wrappers for the orbit class's state variable setters and getters
    @property
    def orbital_freq(self):

        return self.orbit.get_orbital_freq(self.orbit_location)

    @orbital_freq.setter
    def orbital_freq(self, value: np.ndarray):

        self.orbit.set_orbital_freq(self.orbit_location, value, set_by_planet=True)
        self.orbit_change()

    @property
    def eccentricity(self):
        return self.orbit.get_eccentricity(self.orbit_location)

    @eccentricity.setter
    def eccentricity(self, value: np.ndarray):
        self.orbit.set_eccentricity(self.orbit_location, value, set_by_planet=True)
        self.orbit_change()

    @property
    def inclination(self):
        return self.orbit.get_inclination(self.orbit_location)

    @inclination.setter
    def inclination(self, value: np.ndarray):
        self.orbit.set_inclination(self.orbit_location, value, set_by_planet=True)
        self.orbit_change()

    @property
    def semi_major_axis(self):
        return self.orbit.get_orbital_freq(self.orbit_location)

    @semi_major_axis.setter
    def semi_major_axis(self, value: np.ndarray):
        self.orbit.set_semi_major_axis(self.orbit_location, value, set_by_planet=True)
        self.orbit_change()


class BasicWorld(WorldBase):

    type = 'basic'

    def __init__(self, planet_config: dict, automate: bool = True):
        super().__init__(planet_config, automate=automate)

        self.set_geometry(planet_config['radius'], planet_config['mass'])


class Star(BasicWorld):


    type = 'basic'


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

    type = 'burnman'

    def __init__(self, planet_config: dict, burnman_world: burnman.Planet, bm_layers: list, automate: bool = True):

        super().__init__(planet_config, automate=automate)

        # Additional orbit information is needed for potential target worlds
        self.orbit_num = None

        # Store burnman information
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
            layer = construct_layer(layer_name, self, bm_layer, layer_config)
            self.layers_byname[layer_name] = layer
            self.layers.append(layer)
            setattr(self, layer.name.lower(), layer)

        # Constants used in calculations
        self.tidal_susceptibility_inflated = None # = 1.5 * G * M_pri^2 * R_sec^5 -- Set by orbit class

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

        # Dependent variables
        self._tidal_modes = None
        self._tidal_heating_coeffs = None
        self._tidal_ztorque_coeffs = None
        self._global_love = None
        self._tidal_heating = None
        self._tidal_ztorque = None

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

    def orbit_change(self):

        if self.is_spin_sync:
            self._tidal_modes, self._tidal_heating_coeffs, self._tidal_ztorque_coeffs = \
                spin_sync_modes(self.orbital_freq, self.eccentricity, self.inclination)
        else:
            self._tidal_modes, self._tidal_heating_coeffs, self._tidal_ztorque_coeffs = \
                nsr_modes(self.orbital_freq, self.spin_freq, self.eccentricity, self.inclination)

        for layer in self:
            if isinstance(layer, TidalLayer):
                layer.tidal_modes = self._tidal_modes

        self.update_tides()

    def update_tides(self, set_layer_heating: bool = True):
        """ Calculates tidal heating and tidal zTorque based on current state.

        :return: <Tuple[ndarray, ndarray]> Tidal Heating, Tidal zTorque
        """

        layer_loves = [layer.complex_love for layer in self]
        tidal_susceptibility = self.tidal_susceptibility_inflated / self.semi_major_axis**6

        shape = None
        rheos = list()
        for layer_love in layer_loves:
            for rheo_name, mode_list in layer_love.items():
                if rheo_name not in rheos:
                    rheos.append(rheo_name)
                if mode_list[0].shape in [tuple(), (1,)]:
                    # One dimensional array, ignore
                    continue
                else:
                    if shape is None:
                        if shape != mode_list[0].shape:
                            raise BadArrayShape
                    else:
                        shape = mode_list[0].shape

        global_love = {rheo_name: np.zeros(shape, dtype=np.complex) for rheo_name in rheos}
        global_heating = {rheo_name: np.zeros(shape, dtype=np.complex) for rheo_name in rheos}
        global_ztorque = {rheo_name: np.zeros(shape, dtype=np.complex) for rheo_name in rheos}
        for layer, layer_love in zip(self, layer_loves):

            layer_heating = {}
            layer_ztorque = {}

            if isinstance(layer, TidalLayer):
                scale = layer.tidal_scale
            else:
                scale = 1.

            for rheo_name, love_list in layer_love.items():

                layer_heating_bymodes = [mode * (-np.imag(love)) * heating_coeff for mode, love, heating_coeff in
                                         zip(self.tidal_modes, love_list, self.tidal_heating_coeffs)]
                layer_ztorque_bymodes = [(-np.imag(love)) * ztorque_coeff for love, ztorque_coeff in
                                         zip(love_list, self.tidal_ztorque_coeffs)]
                layer_heating_sum = tidal_susceptibility * scale * sum(layer_heating_bymodes)
                layer_ztorque_sum = tidal_susceptibility * scale * sum(layer_ztorque_bymodes)
                layer_heating[rheo_name] = layer_heating
                layer_ztorque[rheo_name] = layer_ztorque_sum

                global_love[rheo_name] += scale * sum(love_list)
                global_heating[rheo_name] += layer_heating_sum
                global_ztorque[rheo_name] += layer_ztorque_sum

            if set_layer_heating:
                layer.tidal_heating = layer_heating

        self._global_love = global_love
        self._tidal_heating = global_heating
        self._tidal_ztorque = global_ztorque


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
    def tidal_modes(self):
        return self._tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating_coeffs(self):
        return self._tidal_heating_coeffs

    @tidal_heating_coeffs.setter
    def tidal_heating_coeffs(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_ztorque_coeffs(self):
        return self._tidal_ztorque_coeffs

    @tidal_ztorque_coeffs.setter
    def tidal_ztorque_coeffs(self, value):
        raise ImproperAttributeHandling

    @property
    def global_love(self):
        return self._global_love

    @global_love.setter
    def global_love(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating(self):
        return self._tidal_heating

    @tidal_heating.setter
    def tidal_heating(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_ztorque(self):
        return self._tidal_ztorque

    @tidal_ztorque.setter
    def tidal_ztorque(self, value):
        raise ImproperAttributeHandling


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