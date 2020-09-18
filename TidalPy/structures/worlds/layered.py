from typing import Tuple, Dict, Union, Iterator

import numpy as np

from .tidal import TidalWorld
from ..layers import PhysicsLayer, layers_class_by_world_class, LayerType, GasLayer, BurnmanLayer
from ... import log
from ...exceptions import ParameterMissingError, TidalPyWorldError, MissingArgumentError, \
    InitiatedPropertyChangeError
from ...utilities.numpyHelper import find_nearest

BAD_LAYER_SYMBOLS = (' ', '*', '-', '/', '+', '=', '@', '#', '$', '%', '\\', '^', '&', '(', ')', '~', '`')

class LayeredWorld(TidalWorld):

    """ LayeredWorld - Construct tidal worlds that have layers.


    See Also
    --------
    Parent Class:
        TidalPy.structures.worlds.TidalWorld
    Child Classes:
        TidalPy.structures.world.GasGiantLayeredWorld
        TidalPy.structures.world.BurnManWorld
    """

    world_class = 'layered'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):
        """ LayeredWorld constructor

        Parameters
        ----------
        world_config : dict
            Configuration file used to build the world. User provided configs override default configurations that
                TidalPy assumes.
            Please see files stored in <TidalPy directory>/structures/worldConfigs for example configuration dict.
        name : str = None
            Name of the world. If None, will use name provided in world_config.
        initialize : bool = True
            Determines if initial reinit should be performed on the world (loading in data from world_config).
        """

        if 'layers' not in world_config:
            log.error("Layered world's configurations do not contain layer information. "
                      "Construction can not be completed.")
            raise ParameterMissingError("Layered world's configurations do not contain layer information.")

        super().__init__(world_config, name, initialize=False)

        # Basic layer reinit
        LayerClass = layers_class_by_world_class[self.world_class]
        self._layers_class = LayerClass.layer_class
        self._num_layers = len(self.config['layers'])
        last_layer_i = self._num_layers - 1

        # Layer storage properties
        self._layers_by_name = dict()
        self._layers = []
        self._layer_types = []

        # Get layer types
        for layer_i, (layer_name, layer_config) in enumerate(self.config['layers'].items()):

            # Check for issues
            if any(bad_symbol in layer_name for bad_symbol in BAD_LAYER_SYMBOLS):
                log.error(f'An illegal symbol was found in {layer_name} for {self}.')
                raise TidalPyWorldError('Illegal symbol found in layer name.')

            # Pull out configuration info
            layer_type = layer_config['type']
            self._layer_types.append(layer_type)

            # Build Layer
            is_last_layer = layer_i == last_layer_i
            layer = LayerClass(layer_name, layer_i, self, layer_config, is_last_layer, initialize=False)

            # Store layer in the world's containers
            self._layers.append(layer)
            self._layers_by_name[layer_name] = layer

            # Also store it in the world itself as its own attribute
            setattr(self, layer_name, layer)

            # Store aliased layer names
            if layer_name != layer_name.lower():
                self._layers_by_name[layer_name.lower()] = layer
                setattr(self, layer_name.lower(), layer)
            if layer_name != layer_name.title():
                self._layers_by_name[layer_name.title()] = layer
                setattr(self, layer_name.title(), layer)

        # Make layer storage immutable
        self._layers = tuple(self._layers)
        self._layer_types = tuple(self._layer_types)

        if initialize:
            self.reinit(initial_init=True, setup_simple_tides=False, reinit_layers=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, setup_simple_tides: bool = False,
               set_by_burnman: bool = False, reinit_layers: bool = True):
        """ Initialize or Reinitialize the world based on changes to its configurations.

        This must be called at least once before an instance can be used. The constructor will automatically make an
            initial call to reinit unless told to not to.

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this function has been called.
        reinit_geometry : bool = True
            If `True`, the initializer will automatically call the `set_geometry()` method.
        set_by_burnman : bool = False
            Set to `True` if called from a burnman world.
        setup_simple_tides : bool = True
            Set to `True` if a global CPL/CTL tidal calculation is desired.
        reinit_layers : bool = True
            If `True`, calls to the world's layers' reinit() method.
        """

        # Don't let parent classes initialize geometry since a LayeredWorld's mass is based on its layers' masses
        super().reinit(initial_init=initial_init, reinit_geometry=False,
                       set_by_burnman=set_by_burnman,
                       setup_simple_tides=setup_simple_tides)

        # Setup Geometry
        if set_by_burnman:
            radius = self.radius
            mass = self.mass
            update_state_geometry = False
        else:
            # Pull out planet configurations
            radius = self.config['radius']
            volume = (4. / 3.) * np.pi * radius**3
            mass = self.config.get('mass', None)
            update_state_geometry = True

            # Layer constructor may need the planets mass and radius.
            #     So set those here (they will be reset by the set_geometry method).
            self._radius = radius
            self._mass = mass
            self._volume = volume

        # Update the global tidal volume fraction
        running_tidal_fraction = 0.
        running_layer_masses = 0.

        # Setup Layers
        if reinit_layers:
            # Call reinit to the layers within this planet building from bottom to top.
            for layer in self.layers:
                layer.reinit(initial_init, initialize_geometry=True)

                running_layer_masses += layer.mass
                if layer.is_tidal:
                    running_tidal_fraction += layer.tidal_scale

            # Record the world's tidal scale as the sum of its layer's scales
            self.tidal_scale = running_tidal_fraction

            # Pull out densities and pressures and convert them into constant value slices
            # Store some layer information at the world-level
            radii = tuple([layer.radii for layer in self])
            volume_slices = tuple([layer.volume_slices for layer in self])
            depths = tuple([layer.depths for layer in self])
            mass_slices = tuple([layer.mass_slices for layer in self])
            mass_below_slices = tuple([layer.mass_below_slices for layer in self])
            density_slices = tuple([layer.density_slices for layer in self])
            gravity_slices = tuple([layer.gravity_slices for layer in self])

            self._radii = np.concatenate(radii)
            self._volume_slices = np.concatenate(volume_slices)
            self._depths = np.concatenate(depths)
            self._mass_slices = np.concatenate(mass_slices)
            self._mass_below_slices = np.concatenate(mass_below_slices)
            self._density_slices = np.concatenate(density_slices)
            self._gravity_slices = np.concatenate(gravity_slices)
            self._num_slices = len(self._radii)

            if self.mass is None:
                mass = running_layer_masses
            else:
                mass = self.mass

            if reinit_geometry:
                self.set_geometry(radius, mass, thickness=None, mass_below=0.,
                                  update_state_geometry=update_state_geometry, build_slices=False)
                reinit_geometry = False

            # Working from the top-most layer downwards, calculate pressures in each layer.
            #    Start with the pressure_above for this world (to account for a high pressure atmo, etc.)
            self.set_static_pressure(self.pressure_above, build_slices=True)

        if reinit_geometry:
            self.set_geometry(radius, mass, thickness=None, mass_below=0., update_state_geometry=update_state_geometry,
                              build_slices=False)

        # Clean up config:
        if 'radii' in self.config:
            del self._config['radii']
        for layer_name, layer_config in self._config['layers'].items():
            if 'radii' in layer_config:
                del layer_config['radii']

    def set_geometry(self, radius: float, mass: float, thickness: float = None, mass_below: float = 0.,
                     update_state_geometry: bool = True, build_slices: bool = False):
        """ Calculates and sets world's physical parameters based on user provided input.

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
        mass_below : float = 0.
            Mass below this object (only applicable for shell-like structures)
            Used in gravity and pressure calculations
        update_state_geometry : bool = True
            Update the class' state geometry
        build_slices : bool = False
            If True, method will attempt to calculate gravities, densities, etc. for each slice.

        """

        super().set_geometry(radius, mass, thickness=None, mass_below=0., update_state_geometry=update_state_geometry,
                             build_slices=False)

        # For a layered world the middle gravity and pressure will depend upon the layer structure and densities.
        #    We must use the layers' static pressure to determine this world's middle pressure
        if self.radii is not None:
            median_radius = float(np.median(self.radii))
            median_index = find_nearest(self.radii, median_radius)
            self._radius_middle = median_radius
            self._gravity_middle = self.gravity_slices[median_index]
            self._density_middle = self.density_slices[median_index]

    def set_static_pressure(self, pressure_above: float = None, build_slices: bool = True,
                            called_from_burnman: bool = False):
        """ Sets the static pressure for the physical structure.

        `Static` here indicates that this is not a dynamic pressure used in many calculations. The static pressure can
            be used in place of the dynamic pressure, but that is not always the case.

        Parameters
        ----------
        pressure_above : float = None
            Pressure above this structure. If this is a layer, then it is the pressure at the base of the overlying
                layer. If it is the upper-most layer or a world, then it may be the surface pressure.
        build_slices : bool = True
            If `True`, method will find the pressure at each slice of the physical object.
        called_from_burnman : bool = False
            Set to `True` if called from a burnman layer/world.
        """

        log.debug(f'Setting up static pressure for {self}.')

        if called_from_burnman:
            log.debug('Static pressure method called from burnman - skipping.')
            return True

        if pressure_above is None:
            pressure_above = self.pressure_above

        if pressure_above is None:
            log.error(f'Not enough information to build static pressure for {self}.')
            raise MissingArgumentError

        else:
            # Calculate pressures from top down
            self._pressure_outer = pressure_above

            # For a layered world the middle pressure will depend upon the layer structure and densities. We must use
            #    the layers' static pressure to determine this world's middle pressure
            self._pressure_middle = None

            if build_slices:
                # Working from the top-most layer downwards, calculate pressures in each layer.
                #    Start with the pressure_above for this world (to account for a high pressure atmo, etc.)
                running_pressure = self.pressure_above
                for layer in self.layers[::-1]:
                    layer.set_static_pressure(pressure_above=running_pressure, build_slices=True)
                    running_pressure = layer.pressure_inner

                # Now pull out pressure slice data
                pressure_slices = tuple([layer.pressure_slices for layer in self])
                self._pressure_slices = np.concatenate(pressure_slices)

                # Find the median radius and use it to find the "middle" pressure
                median_radius = float(np.median(self.radii))
                median_index = find_nearest(self.radii, median_radius)
                self._pressure_middle = self.pressure_slices[median_index]
                self._pressure_inner = self.pressure_slices[0]

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

    def update_time(self):
        """ Update any attributes or methods that rely on the current time.

        This can be called by a Orbit class if one is set.
        """

        super().update_time()

        # Tell each layer that the time has been updated.
        for layer in self:
            layer.update_time()

    def update_orbit_spin(self, frequency_change: bool = True):
        """ Performs state updates whenever the planet's orbital parameters are changed

        Parameters
        ----------
        frequency_change : bool = True
            If the orbit or spin rate changed then additional calculations must take place.
        """

        if frequency_change:
            # A change to the orbital or spin frequency will cause a change to the forcing frequency(ies) which will
            #    change the complex compliances
            for layer in self:
                layer.rheology.update_frequency(called_from_world=True)

        super().update_orbit_spin(frequency_change=frequency_change)

    # # Initiated properties
    @property
    def layers_class(self) -> str:
        """ Name of the python class used to construct this world's layers """
        return self._layers_class

    @layers_class.setter
    def layers_class(self, value):
        raise InitiatedPropertyChangeError

    @property
    def layers_by_name(self) -> Dict[str, LayerType]:
        """ Dictionary of layers with the keys equaling the layer names """
        return self._layers_by_name

    @layers_by_name.setter
    def layers_by_name(self, value):
        raise InitiatedPropertyChangeError

    @property
    def layers(self) -> Tuple[LayerType, ...]:
        """ Tuple of layers within world, ordered from bottom-most to top-most """
        return self._layers

    @layers.setter
    def layers(self, value):
        raise InitiatedPropertyChangeError

    @property
    def layer_types(self) -> Tuple[str, ...]:
        """ Tuple of layer types within world, ordered from bottom-most to top-most """
        return self._layer_types

    @layer_types.setter
    def layer_types(self, value):
        raise InitiatedPropertyChangeError

    @property
    def num_layers(self) -> int:
        """ Number of layers within world """
        return self._num_layers

    @num_layers.setter
    def num_layers(self, value):
        raise InitiatedPropertyChangeError


    # # Aliased properties
    @property
    def densities(self):
        return self.density_slices

    @densities.setter
    def densities(self, value):
        self.density_slices = value

    @property
    def gravities(self):
        return self.gravity_slices

    @gravities.setter
    def gravities(self, value):
        self.gravity_slices = value

    @property
    def pressures(self):
        return self.pressure_slices

    @pressures.setter
    def pressures(self, value):
        self.pressure_slices = value


    # Dunder properties
    def __iter__(self) -> Iterator[Union[PhysicsLayer, GasLayer, BurnmanLayer]]:
        """ Planet will iterate over its layers
        Returns
        -------
        iter(self.layers) : Iterator[Union[PhysicsLayer, GasLayer, BurnmanLayer]]
            The iterator of the layer list.
        """

        return iter(self.layers)
