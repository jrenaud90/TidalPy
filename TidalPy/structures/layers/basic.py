from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

from .defaults import layer_defaults
from .helper import find_geometry_from_config
from ..physical import PhysicalObjSpherical
from ... import log
from ...exceptions import OuterscopePropertySetError, MissingArgumentError, \
    ConfigPropertyChangeError, ParameterMissingError, InitiatedPropertyChangeError, IncorrectMethodToSetStateProperty
from ...utilities.types import NoneType, FloatArray

if TYPE_CHECKING:
    from ..worlds import LayeredWorldType
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
    TidalPy.structures.layers.BurnmanLayer
    TidalPy.structures.layers.GasLayer
    """

    default_config = layer_defaults
    layer_class = 'base'

    def __init__(self, layer_name: str, layer_index: int, world: 'LayeredWorldType', layer_config: dict,
                 is_top_layer: bool, initialize: bool = True):
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
        self.default_config = self.default_config[self.type]

        # Initiated Attributes
        self._name = layer_name
        self._layer_index = layer_index
        self._world = world
        self._is_top_layer = is_top_layer

        # Setup Physical Layer Geometry
        super().__init__(layer_config)

        # State properties
        self._gravity = None
        self._density = None
        self._pressure = None
        self._temperature = None

        # Other attributes
        self.material_name = self.config['material']
        self.tidal_scale = 1.
        self.heat_sources = None

        # Flags
        self._is_tidal = None
        self._use_tidal_vol_frac = None

        log.debug(f'Creating: {self}.')
        if initialize:
            self.reinit(initial_init=initialize)

    def reinit(self, initial_init: bool = False, set_by_burnman: bool = False, initialize_geometry: bool = True):
        """ Reinitialize the physical object by pulling in any potentially new configurations

        Parameters
        ----------
        initial_init : bool = False
            Set to `True` for the first time an instance is created.
        set_by_burnman : bool = False
            Set to `True` if a Burnman layer/world constructor is calling reinit
        initialize_geometry : bool = False
            Set to `True` if the set_geometry method should be called from within reinit
        """

        super().reinit(initial_init=initial_init, set_by_burnman=set_by_burnman)

        # Load in configurations
        self._is_tidal = self.config['is_tidally_active']
        self._use_tidal_vol_frac = self.config['use_tidal_vol_frac']

        if self.config['use_surf_gravity']:
            # Use surface gravity for layer instead of the gravity set by interpolating burnman data (mid/avg/etc)
            # This primarily affects convection calculation
            self._gravity = lambda: self.gravity_surf

        if not set_by_burnman and initialize_geometry:
            try:
                radius_layer_below = None
                if self.layer_below is not None:
                    radius_layer_below = self.layer_below.radius

                radius, thickness, volume, mass, density = \
                    find_geometry_from_config(self.config, self.layer_index, self.is_top_layer,
                                              self.world.radius, self.world.mass, radius_layer_below)
            except ParameterMissingError as e:
                log.error(f'Not enough information provided to determine geometry for {self}.')
                raise e

            # Set the layer's geometry
            #     OPT: some of set_geometry will end up redoing some of the above calculations, but they are not
            #         expensive calculations and it should only be preformed a handful of times.
            #         But perhaps an area for future optimization.
            self.set_geometry(radius=radius, mass=mass, thickness=thickness, update_state_geometry=True,
                              build_slices=True)

        # Clean up config:
        if 'radii' in self.config:
            del self._config['radii']

    def set_geometry(self, radius: float, mass: float, thickness: float = None,
                     mass_below: float = None, update_state_geometry: bool = True, build_slices: bool = True):
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
                log.error(f'Layer thickness can not be determined for {self}')
                raise MissingArgumentError

        # Gravity (and therefore pressure) will depend on the mass contained below this layer
        if self.layer_index == 0 :
            # Bottom-most layer - nothing below it.
            mass_below = 0.
        else:
            mass_below = sum([self.world.layers[i].mass for i in range(0, self.layer_index)])

        # Setup Layer Geometry
        super().set_geometry(radius, mass, thickness, mass_below=mass_below,
                             update_state_geometry=update_state_geometry, build_slices=build_slices)

        # Setup Tidal Volume Fraction
        if self.use_tidal_vol_frac:
            self.tidal_scale = self.volume / self.world.volume

    def update_thermal(self):
        """ Update various classes and methods when any thermal parameters change. """

        log.debug(f'Update thermals called for {self}.')

        # Thermal updates are done by child classes

    def update_time(self):
        """ Update various classes and methods when the time of the world changes. """

        log.debug(f'Update time called for {self}')

        # Time updates are done by child classes

    def clear_state(self, clear_pressure: bool = False):

        log.debug(f'Clear state called for {self}. Clear pressure = {clear_pressure}.')

        super().clear_state()

        self._temperature = None
        if clear_pressure:
            self._pressure = None

    def set_state(self, temperature: FloatArray = None, pressure: FloatArray = None, call_updates: bool = True):
        """ Set the layer's state properties

        Parameters
        ----------
        temperature : FloatArray
            New dynamic temperature for the layer.
        pressure :
            New dynamic pressure for the layer.
        call_updates : bool = True
            If `True`, method will call the update thermals method.

        """

        # Check if temperature or pressure were provided, call the respective methods but hold off on updating until
        #    the end (increase to performance).

        if temperature is not None:
            self.set_temperature(temperature, call_updates=False)

        if pressure is not None:
            self.set_pressure(pressure, call_updates=False)

        if call_updates:
            self.update_thermal()

    def set_temperature(self, temperature : FloatArray, call_updates: bool = True):
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
            self.update_thermal()

    def set_pressure(self, pressure: FloatArray, call_updates: bool = True):
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
            self.update_thermal()


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

    # TODO: How do we want to handle these?
    @property
    def gravity(self) -> float:
        """ State gravity of the layer """
        return self._gravity

    @gravity.setter
    def gravity(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def density(self) -> float:
        """ State density of the layer """
        return self._density

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
    def temperature_surf(self):
        """ The temperature at the top of this layer """
        if self.is_top_layer:
            # Return the world's surface temperature
            return self.world.surface_temperature
        else:
            # Return the dynamic temperature of the layer above it
            return self.layer_above.temperature

    @temperature_surf.setter
    def temperature_surf(self, value):
        raise OuterscopePropertySetError


    # # Aliased properties
    @property
    def surface_temperature(self):
        """ Alias for self.temperature_surf """
        return self.temperature_surf

    @surface_temperature.setter
    def surface_temperature(self, value):
        self.temperature_surf = value


    # # Dunder methods
    def __str__(self):

        if self.world is None:
            text = f'[Layer {self.name} ({self.type}) no world]'
        else:
            text = f'[Layer {self.name} ({self.type}) in {self.world} (loc={self.layer_index})]'
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
