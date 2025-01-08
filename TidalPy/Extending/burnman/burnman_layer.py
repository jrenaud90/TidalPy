from typing import TYPE_CHECKING

burnman_installed = True
try:
    import burnman
except ImportError:
    burnman_installed = False
    # Build fake class so type checking passes.
    class burnman:
        Planet = None
        Layer = None
        Material = None

import numpy as np


import TidalPy
from TidalPy.exceptions import (AttributeNotSetError, ImproperPropertyHandling, UnknownTidalPyConfigValue)
from TidalPy.utilities.numpy_helper.array_other import find_nearest
from TidalPy.utilities.types import FloatArray
from TidalPy.structures.layers.physics import PhysicsLayer
from TidalPy.Extending.burnman.conversion import burnman_property_name_conversion, burnman_property_value_conversion

from TidalPy.logger import get_logger
log = get_logger("TidalPy")

if TYPE_CHECKING:
    from .burnman_world import BurnManWorld


class BurnmanLayer(PhysicsLayer):
    """ BurnmanLayer
    Layer object to store parameters calculated by the Burnman software. Additionally it contains properties and
        functionality that matches the TidalPy `PhysicsLayer`.


    See Also
    --------
    TidalPy.structures.layers.PhysicsLayer
    """

    layer_class = 'burnman'

    def __init__(
        self, layer_name: str, layer_index: int, world: 'BurnManWorld', layer_config: dict,
        is_top_layer: bool, initialize: bool = True
        ):
        """ Burnman layer constructor

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

        if not burnman_installed:
            log.warn('Burnman package not found. BurnmanLayer will have limited functionality.')

        super().__init__(layer_name, layer_index, world, layer_config, is_top_layer, initialize=False)

        # BurnmanLayer is built on the BurnMan API and layer system.
        #    Pull out BurnMan parameters
        burnman_layer = self.world.bm_layers[layer_index]
        self._bm_layer = burnman_layer
        self._bm_material = burnman_layer.material

        # Other state variables that will be initialized later
        self._bm_mid_index = None
        self._persistent_pressure = None
        self.interp_func = None

        # Material Properties set by BurnMan
        self.interp_temperature_range = np.linspace(
            *tuple(self.config['interp_temperature_range']),
            self.world.config['bm_interpolation_n']
            )
        self._interp_prop_data_lookup = dict()

        if initialize:
            self.reinit(initial_init=initialize)

    def reinit(self, initial_init: bool = False, initialize_geometry: bool = False):
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

        if initial_init or initialize_geometry:
            # Find how many slices
            self._num_slices = len(self.bm_layer.radii)

            # Find which slice index is nearest to the middle of the layer (as determined by the radius/thickness)
            thickness = self.bm_layer.thickness
            self._bm_mid_index = find_nearest(self.bm_layer.radii, self.bm_layer.radii - (thickness / 2.))

            # Setup geometry
            self.set_geometry(radius=self.bm_layer.outer_radius, mass=self.bm_layer.mass, thickness=thickness)

        # Call the standard reinit for the physics layer
        super().reinit(initial_init, initialize_geometry=False)

        # Build lookup tables
        self._build_material_property_interpolation()

    def set_state(
        self, temperature: FloatArray = None, pressure: FloatArray = None, viscosity: FloatArray = None,
        shear_modulus: FloatArray = None
        ):
        """ Set the layer's state properties

        Parameters
        ----------
        temperature : FloatArray = None
            New dynamic temperature for the layer [K].
        pressure : FloatArray = None
            New dynamic pressure for the layer [Pa].
        viscosity : FloatArray = None
            The new viscosity of the layer in [Pa s]
        shear_modulus : FloatArray = None
            The new shear modulus of the layer in [Pa]

        """

        # Check if the pressure or temperature of the layer has changed.
        temperature_change = temperature is not None
        pressure_change = pressure is not None

        # TODO: pressure_change is currently unused (see TODO below).

        if temperature_change:
            # Update material properties interpolated by burnman.

            # TODO: Current limitations:
            #    - No material properties change when the pressure changes.
            #        - Future: Build another look up table for pressure?
            #    - Does not make calls to Burnman's EOS calculator for temperature changes, instead a lookup table is
            #       used to find pre-calculated values for the material properties.
            #        - Generally this should be okay for small temperature changes, but the larger the delta-temp the
            #           larger the error will be.
            #    - Recalculating the EOS can be computationally expensive, especially once it is wrapped within a bunch
            #       of TidalPy methods. This may be a major limitation to TidalPy's future development.
            #    - Workaround: the user can access the burnman layer and planet to manually update the EOS, but changes
            #       won't propagate into the TidalPy methods (as of writing).

            # Set new material properties based on BurnMan Interpolation
            for prop_name, prop_data in self._interp_prop_data_lookup.items():
                setattr(self, prop_name, np.interp(temperature, self.interp_temperature_range, prop_data))
            self.thermal_diffusivity = self.thermal_conductivity / (self.density * self.specific_heat)

        # Set the new state properties in the layer.
        super().set_state(temperature=temperature, pressure=pressure, viscosity=viscosity, shear_modulus=shear_modulus)

    def set_geometry(
        self, radius: float = None, mass: float = None, thickness: float = None,
        mass_below: float = 0., update_state_geometry: bool = True, build_slices: bool = True
        ):
        """ Calculates and sets the layer's physical parameters based on user provided input.

        For a BurnmanLayer almost all of the geometry is set by the results of the burnman evaluation.

        Assumptions
        -----------
        Spherical Geometry

        Parameters
        ----------
        radius : float = None
            Outer radius of object [m]
        mass : float = None
            Mass of object [kg]
        thickness : float = None
            Thickness of the object [m]
        mass_below : float = 0.
            Mass below this object (only applicable for shell-like structures)
            Used in gravity and pressure calculations
        update_state_geometry : bool = True
            Update the class' state geometry
        build_slices : bool = True
            If True, method will attempt to calculate gravities, densities, etc. for each slice.

        """

        # Set the physical properties
        self._radius = self.bm_layer.outer_radius
        self._thickness = self.bm_layer.thickness
        self._radius_middle = self.radius - (self.thickness / 2.)
        self._radius_inner = self.radius - self.thickness
        self._volume = (4. / 3.) * np.pi * (self.radius**3 - self.radius_inner**3)
        self._surface_area_outer = 4. * np.pi * self.radius**2
        self._surface_area_middle = 4. * np.pi * self.radius_middle**2
        self._surface_area_inner = 4. * np.pi * self.radius_inner**2

        # Set the mass properties
        self._mass = self.bm_layer.mass
        self._density_bulk = self.mass / self.volume
        if self.layer_below is None:
            layer_below_mass = 0.
        else:
            layer_below_mass = self.layer_below.mass
        self._mass_below = layer_below_mass

        # Set slice physical properties
        #    For Burnman layers we want to skip the first step, we will do this for all items pulled from BM.
        #    But we want to keep the same number of shells. So we will perform an interpolation to expand from the
        #    Burnman values
        bm_radii = self.bm_layer.radii
        self._radii = np.linspace(bm_radii[1], bm_radii[-1], self.num_slices)
        self._depths = self.radius - self._radii
        self._sa_slices = 4. * np.pi * self.radii**2
        self._volume_slices = np.zeros_like(self._radii)
        self._volume_slices[0] = (4. / 3.) * np.pi * (self.radii[0]**3 - self.radius_inner**3)
        self._volume_slices[1:] = (4. / 3.) * np.pi * (self.radii[1:]**3 - self.radii[:-1]**3)

        # Set slice mass properties
        self._density_slices = np.interp(self.radii, bm_radii, self.bm_layer.density)
        self._gravity_slices = np.interp(self.radii, bm_radii, self.bm_layer.gravity)
        self._pressure_slices = np.interp(self.radii, bm_radii, self.bm_layer.pressures)
        self._mass_slices = self.density_slices * self.volume_slices

        # Mass below each slice is equal to slice masses + and mass below this physical object
        self._mass_below_slices = np.asarray(
            [self.mass_below + sum(self.mass_slices[:i + 1]) for i in range(self.num_slices)]
            )

        # For BurnmanLayer, assume everything is calculated by BurnMan.
        #     Override some of the items set by the set_geometry method.
        if self.world.config['bm_interpolation_method'] == 'mid':
            self._pressure_middle = self.bm_layer.pressures[self._bm_mid_index]
            self._density_middle = self.bm_layer.density[self._bm_mid_index]
            self._gravity_middle = self.bm_layer.gravity[self._bm_mid_index]
            self.interp_func = lambda array: array[self._bm_mid_index]
        elif self.world.config['bm_interpolation_method'] == 'avg':
            self._pressure_middle = np.average(self.bm_layer.pressures)
            self._density_middle = np.average(self.bm_layer.density)
            self._gravity_middle = np.average(self.bm_layer.gravity)
            self.interp_func = np.average
        elif self.world.config['bm_interpolation_method'] == 'median':
            self._pressure_middle = np.median(self.bm_layer.pressures)
            self._density_middle = np.median(self.bm_layer.density)
            self._gravity_middle = np.median(self.bm_layer.gravity)
            self.interp_func = np.median
        else:
            raise UnknownTidalPyConfigValue

        # Set up a pressure that will persist if the layer's state is cleared
        self._persistent_pressure = np.copy(self.pressure)

        # Update state properties
        #    Note: these are not actually accessible without using the protected getter.
        #    Instead, the regular getter method pulls from the BurnMan layer directly (they have been inner-scoped)
        #    These definitions below are set just so that the protected properties still hold a value.
        self._temperature_outer = self.bm_layer.temperatures[-1]
        self._temperature_inner = self.bm_layer.temperatures[0]
        self._density_outer = self.bm_layer.density[-1]
        self._density_inner = self.bm_layer.density[0]
        self._gravity_outer = self.bm_layer.gravity[-1]
        self._gravity_inner = self.bm_layer.gravity[0]
        self._pressure_outer = self.bm_layer.pressures[-1]
        self._pressure_inner = self.bm_layer.pressures[0]

        # Make a call to the parent methods set_geometry method, but update_state_geometry=False will prevent them from
        #    changing the geometry that was set above.
        super().set_geometry(
            self.radius, self.mass, self.thickness, mass_below=layer_below_mass,
            update_state_geometry=False, build_slices=False
            )

    def set_static_pressure(
        self, pressure_above: float = None, build_slices: bool = True,
        called_from_burnman: bool = True
        ):
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
        log.debug('Static pressure method called from burnman layer - skipping update.')
        # super().set_static_pressure(pressure_above, build_slices)

    def geotherm(self, avg_temperature: float = None):
        """ Calculates layer's geotherm based on an average temperature (or layer's current temperature)

        Returns
        -------
        temperature_profile : np.ndarray
            numpy array of the adiabatic temperature profile
        """

        if avg_temperature is None:
            avg_temperature = self.temperature

        temperature_profile = \
            burnman.geotherm.adiabatic(self.bm_layer.pressures, avg_temperature, self.bm_layer.material)

        # We again want to shrink the temperature domain to match the radii which only take [1:] of the Burnman version.
        bm_radii = self.bm_layer.radii
        temperature_profile = np.interp(self.radii, bm_radii, temperature_profile)

        return temperature_profile

    def _build_material_property_interpolation(self):
        """ Interpolates material properties based on a fixed pressure and a suggested temperature range.

        Using the BurnMan package's equation of states, this function will build interpolated lookup tables for several
            material properties. These lookup tables are functions of only temperature; pressure is assumed to not
            change significantly in a TidalPy simulation. However, the machinery to change properties based on
            pressure is built into BurnMan and can be accessed via the layer's reference to the Burnman layer:
                self.bm_material.evaluate

        The specific constant pressure used in building the lookup table is set by the layer configuration. See the
            variable 'bm_interpolation_method' in TidalPy's main configurations.py

        Conversions between BurnMan and TidalPy property names are also performed here. Some conversions require a
             transition from molar to specific.

        The final interpolated lookup table is stored in the self._interp_prop_data_lookup which is used in the
             temperature.setter
        """

        log.debug(f'Building Burnman material property lookup tables for {self}.')

        if self.pressure is None:
            raise AttributeNotSetError

        self._interp_prop_data_lookup = dict()
        interp_properties = ['bulk_modulus', 'thermal_expansion', 'specific_heat']
        bm_properties = [burnman_property_name_conversion[interp_prop] for interp_prop in interp_properties]

        # We will use the fixed pressure to calculate the various parameters at all the temperature ranges
        #     These results will then be used in an interpolation for whenever self.temperature changes
        pressures = self.pressure * np.ones_like(self.interp_temperature_range)
        property_results = self.bm_material.evaluate(bm_properties, pressures, self.interp_temperature_range)

        for interp_prop, prop_result in zip(interp_properties, property_results):

            # Perform any unit or other conversions needed to interface BurnMan to TidalPy
            conversion_type = burnman_property_value_conversion.get(interp_prop, None)
            if conversion_type is not None:
                if conversion_type == 'molar':
                    # Convert the parameter from a molar value to a specific value
                    prop_result = prop_result / self.bm_material.molar_mass
                else:
                    raise KeyError

            self._interp_prop_data_lookup[interp_prop] = np.asarray(prop_result)

    # State properties
    @property
    def bm_layer(self) -> burnman.Layer:
        return self._bm_layer

    @bm_layer.setter
    def bm_layer(self, value):
        raise ImproperPropertyHandling('Can not change burnman layer information after layer has been initialized.')

    @property
    def bm_material(self) -> burnman.Material:
        return self._bm_material

    @bm_material.setter
    def bm_material(self, value):
        raise ImproperPropertyHandling('Can not change burnman layer information after layer has been initialized.')

    @property
    def pressure(self):
        """ Layer's dynamic pressure used in calculations """
        # If there is no user provided pressure, use the one set by burnman
        if self._pressure is None:
            return self.pressure_middle
        else:
            return self._pressure

    @pressure.setter
    def pressure(self, new_pressure: FloatArray):
        self.set_pressure(new_pressure)

# Add burnman layer type to known layer types
log.debug('Adding BurnmanLayer to known layer classes.')
TidalPy.structures.layers.known_layer_classes['burnman'] = BurnmanLayer
log.debug('Adding BurnmanLayer to known layer classes (by world type).')
TidalPy.structures.layers.layers_class_by_world_class['burnman'] = BurnmanLayer
