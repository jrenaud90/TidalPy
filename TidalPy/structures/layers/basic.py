from __future__ import annotations

from typing import TYPE_CHECKING, Union

import TidalPy
from TidalPy.exceptions import (ConfigPropertyChangeError, IncorrectMethodToSetStateProperty,
                                InitiatedPropertyChangeError, MissingArgumentError, OuterscopePropertySetError,
                                ParameterMissingError)
from .helper import find_geometry_from_config
from ..physical import PhysicalObjSpherical

from TidalPy.logger import get_logger
log = get_logger("TidalPy")

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray, NoneType

    from ..world_types import LayeredWorldType
    from . import LayerType


class LayerBase(PhysicalObjSpherical):
    """ LayerBase
    Layer object to store parameters geometric and physical properties calculated by TidalPy based on a user-provided
        configuration dictionary.

    Notes:
    .. Does not provide any functionality to perform tidal calculations (see PhysicsLayer instead)

    See Also
    --------
    TidalPy.structures.layers.PhysicsLayer
    TidalPy.structures.layers.GasLayer
    """

    layer_class = 'base'

    def __init__(
        self, layer_name: str, layer_index: int, world: 'LayeredWorldType', layer_config: dict,
        is_top_layer: bool, initialize: bool = True
        ):
        """ Basic layer constructor

        Parameters
        ----------
        layer_name : str
            User-friendly name of layer.
        layer_index : int
            Location of layer within a world (0 indicates center-most).
        world : LayeredWorldType
            World instance where layer was initialized in.
        layer_config : dict
            Layer's user-provided configurations.
        is_top_layer : bool
            If `True`, this layer is the top-most layer.
        initialize : bool = True
            If `True`, then the Layer's reinit is called at the end of the constructor.
        """

        # Load layer defaults based on layer type
        self.type = layer_config['type']
        self.default_config = TidalPy.config['layers'][self.type]

        # Initiated Attributes
        self._name = layer_name
        self._layer_index = layer_index
        self._world = world
        self._is_top_layer = is_top_layer

        # Setup Physical Layer Geometry
        super().__init__(layer_config)

        # State properties
        self._pressure = None
        self._temperature = None

        # Other attributes
        self.material_name = self.config['material']
        self.tidal_scale = 1.
        self.heat_sources = None

        # Configuration properties
        self._is_tidal = None
        self._use_tidal_vol_frac = None
        self._use_surf_gravity = None
        self._use_bulk_density = None

        log.debug(f'Creating: {self}.')
        if initialize:
            self.reinit(initial_init=initialize)

    def reinit(self, initial_init: bool = False, initialize_geometry: bool = True):
        """ Reinitialize the physical object by pulling in any potentially new configurations

        Parameters
        ----------
        initial_init : bool = False
            Set to `True` for the first time an instance is created.
        initialize_geometry : bool = False
            Set to `True` if the set_geometry method should be called from within reinit
        """

        super().reinit(initial_init=initial_init)

        # Load in configurations
        self._is_tidal = self.config['is_tidally_active']
        self._use_tidal_vol_frac = self.config['use_tidal_vol_frac']
        self._use_surf_gravity = self.config['use_surface_gravity']
        self._use_bulk_density = self.config['use_bulk_density']

        if initialize_geometry:
            radius_layer_below = None
            if self.layer_below is not None:
                radius_layer_below = self.layer_below.radius

            radius, thickness, volume, mass, density = \
                find_geometry_from_config(
                    self.config, self.layer_index, self.is_top_layer,
                    self.world.radius, self.world.mass, radius_layer_below
                    )

            # Set the layer's geometry
            #     OPT: some of set_geometry will end up redoing some of the above calculations, but they are not
            #         expensive calculations and it should only be preformed a handful of times.
            #         But perhaps an area for future optimization.
            self.set_geometry(
                radius=radius, mass=mass, thickness=thickness, update_state_geometry=True,
                build_slices=True
                )

        # Clean up config:
        if 'radii' in self.config:
            del self._config['radii']

    def time_changed(self):
        """ The time has changed. Make any necessary updates. """

        log.debug(f'Time changed called for {self}.')

        # Updated set by child methods

    def internal_thermal_equilibrium_changed(self):
        """ The internal heating / cooling of the layer has changed. Make any necessary updates. """

        log.debug(f'Internal thermal equilibrium change called for {self}.')

    def surface_temperature_changed(self, called_from_cooling: bool = False):
        """ Surface temperature has changed - Perform any calculations that may have also changed.

        Parameters
        ----------
        called_from_cooling : bool = False
            Flag to avoid recursive loops between surface temperature and cooling.
        """

        log.debug(f'Surface temperature changed called for {self}.')

        # Updates set by child methods.

    def tidal_frequencies_changed(self, collapse_tidal_modes: bool = True):
        """ The tidal frequencies have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        log.debug(f'Tidal frequencies changed for {self}.')

        # Updates set by child methods.

    def temperature_pressure_changed(self):
        """ The temperature and/or pressure of the layer has changed. Make any necessary updates. """

        log.debug(f'Temperature and/or pressure changed for {self}.')

        # Updates set by child methods.

    def strength_changed(self):
        """ The viscosity and/or shear modulus of the layer has changed. Make any necessary updates. """

        log.debug(f'Strength changed for {self}.')

        # Updates set by child methods.

    def clear_state(self, clear_pressure: bool = False):

        log.debug(f'Clear state called for {self}. Clear pressure = {clear_pressure}.')

        super().clear_state()

        self._temperature = None
        if clear_pressure:
            self._pressure = None

    def set_state(self, temperature: 'FloatArray' = None, pressure: 'FloatArray' = None):
        """ Set the layer's state properties

        Parameters
        ----------
        temperature : FloatArray = None
            New dynamic temperature for the layer [K].
        pressure : FloatArray = None
            New dynamic pressure for the layer [Pa].

        """

        # Check if temperature or pressure were provided, call the respective methods but hold off on updating until
        #    the end (increase to performance).
        temp_or_press_changed = False

        if temperature is not None:
            self.set_temperature(temperature, call_updates=False)
            temp_or_press_changed = True

        if pressure is not None:
            self.set_pressure(pressure, call_updates=False)
            temp_or_press_changed = True

        if temp_or_press_changed:
            self.temperature_pressure_changed()

    def set_geometry(
        self, radius: float, mass: float, thickness: float = None,
        mass_below: float = None, update_state_geometry: bool = True, build_slices: bool = True
        ):
        """ Calculates and sets the layer's physical parameters based on user provided input.

        Assumptions
        -----------
        Spherical Geometry

        Parameters
        ----------
        radius : float
            Outer radius of object [m]
        mass : float
            Mass of object [kg]
        thickness : float = None
            Thickness of the object [m]
        mass_below : float = None
            NOT USED: Left here to keep the method's signature the same as its parents.
            The method will determine for itself what the mass below is based on `self.layer_below` which is a wrapper
                to the layer's world class.
            Mass below this object (only applicable for shell-like structures)
            Used in gravity and pressure calculations
        update_state_geometry : bool = True
            Update the class' state geometry
        build_slices : bool = False
            If True, method will attempt to calculate gravities, densities, etc. for each slice.

        """

        if thickness is None:
            if self.layer_index == 0:
                # Bottom layer: thickness = radius
                thickness = radius
            else:
                raise MissingArgumentError

        # Gravity (and therefore pressure) will depend on the mass contained below this layer
        if self.layer_index == 0:
            # Bottom-most layer - nothing below it.
            mass_below = 0.
        else:
            mass_below = sum([self.world.layers[i].mass for i in range(0, self.layer_index)])

        # Setup Layer Geometry
        super().set_geometry(
            radius, mass, thickness, mass_below=mass_below,
            update_state_geometry=update_state_geometry, build_slices=build_slices
            )

        # Setup Tidal Volume Fraction
        if self.use_tidal_vol_frac:
            self.tidal_scale = self.volume / self.world.volume

    def set_temperature(self, temperature: 'FloatArray', call_updates: bool = True):
        """ Set the layer's dynamic temperature

        Parameters
        ----------
        temperature : FloatArray
            New dynamic temperature for the layer.
        call_updates : bool = True
            If `True`, method will call the update thermals method.

        """

        self._temperature = temperature

        if call_updates:
            self.temperature_pressure_changed()

    def set_pressure(self, pressure: 'FloatArray', call_updates: bool = True):
        """ Set the layer's dynamic pressure

        Parameters
        ----------
        pressure : FloatArray
            New dynamic pressure for the layer.
        call_updates : bool = True
            If `True`, method will call the update thermals method.
        """

        self._pressure = pressure

        if call_updates:
            self.temperature_pressure_changed()

    # # Initialized properties
    @property
    def name(self) -> str:
        """ Name of the layer """
        return self._name

    @name.setter
    def name(self, value):
        raise InitiatedPropertyChangeError

    @property
    def layer_index(self) -> int:
        """ Index of where the layer is inside of `<layer>.world` """
        return self._layer_index

    @layer_index.setter
    def layer_index(self, value):
        raise InitiatedPropertyChangeError

    @property
    def world(self) -> 'LayeredWorldType':
        """ The world class where this layer was initialized """
        return self._world

    @world.setter
    def world(self, value):
        raise InitiatedPropertyChangeError

    @property
    def is_top_layer(self) -> bool:
        """ Flag for if this layer is the top-most layer inside of `<layer>.world` """
        return self._is_top_layer

    @is_top_layer.setter
    def is_top_layer(self, value):
        raise InitiatedPropertyChangeError

    # # State properties
    @property
    def temperature(self) -> FloatArray:
        """ Dynamic layer temperature (taken to be at the interpolation point) [K] """
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self.set_temperature(value)

    @property
    def pressure(self) -> FloatArray:
        """ Dynamic layer pressure (taken to be at the interpolation point) [Pa] """
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        self.set_pressure(value)

    @property
    def gravity(self) -> float:
        """ State gravity of the layer """
        if self._use_surf_gravity:
            return self.gravity_outer
        else:
            return self.gravity_middle

    @gravity.setter
    def gravity(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def density(self) -> float:
        """ State density of the layer """
        if self.use_bulk_density:
            return self.density_bulk
        else:
            return self.density_middle

    @density.setter
    def density(self, value):
        raise IncorrectMethodToSetStateProperty

    # # Configuration Properties
    @property
    def is_tidal(self) -> bool:
        """ Flag if layer is tidally active """
        return self._is_tidal

    @is_tidal.setter
    def is_tidal(self, value):
        raise ConfigPropertyChangeError

    @property
    def use_tidal_vol_frac(self) -> bool:
        """ Flag if layer uses tidal volume fraction """
        return self._use_tidal_vol_frac

    @use_tidal_vol_frac.setter
    def use_tidal_vol_frac(self, value):
        raise ConfigPropertyChangeError

    @property
    def use_surf_gravity(self) -> bool:
        """ Flag if layer uses the surface or central gravity for calculations """
        return self._use_surf_gravity

    @use_surf_gravity.setter
    def use_surf_gravity(self, value):
        raise ConfigPropertyChangeError

    @property
    def use_bulk_density(self) -> bool:
        """ Flag if layer uses the bulk or central density for calculations (only matters for burnman layers) """
        return self._use_bulk_density

    @use_bulk_density.setter
    def use_bulk_density(self, value):
        raise ConfigPropertyChangeError

    # # Outer-scope properties
    #    World Class
    @property
    def layer_below(self) -> Union['LayerType', NoneType]:
        """ The layer instance below this one in a world """
        if self.layer_index == 0:
            return None
        else:
            return self.world.layers[self.layer_index - 1]

    @layer_below.setter
    def layer_below(self, value):
        raise OuterscopePropertySetError

    @property
    def layer_above(self) -> Union['LayerType', NoneType]:
        """ The layer instance above this one in a world """
        if self.layer_index == self.world.num_layers - 1:
            return None
        else:
            return self.world.layers[self.layer_index + 1]

    @layer_above.setter
    def layer_above(self, value):
        raise OuterscopePropertySetError

    @property
    def time(self):
        """ Time property used for radiogenic calculations [Myr]

        Stored in <world>.time
        """

        return self.world.time

    @time.setter
    def time(self, value):
        raise OuterscopePropertySetError

    @property
    def surface_temperature(self):
        """ The temperature at the top of this layer """
        if self.is_top_layer:
            # Return the world's surface temperature
            return self.world.surface_temperature
        else:
            # Return the dynamic temperature of the layer above it
            return self.layer_above.temperature

    @surface_temperature.setter
    def surface_temperature(self, value):
        raise OuterscopePropertySetError

    # # Dunder methods
    def __str__(self):

        if self.world is None:
            text = f'Layer {self.name} ({self.type} no world)'
        else:
            text = f'Layer {self.name} ({self.type} in {self.world}; loc={self.layer_index})'
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
