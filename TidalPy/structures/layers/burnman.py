from typing import TYPE_CHECKING

import burnman
import numpy as np

from .physics import PhysicsLayer
from ...burnman_interface.conversion import burnman_property_name_conversion, burnman_property_value_conversion
from ...configurations import configurations
from ...exceptions import (UnknownTidalPyConfigValue, AttributeNotSetError, ImproperPropertyHandling,
                           IncorrectAttributeType)
from ...utilities.numpyHelper.array_other import find_nearest
from ...utilities.types import FloatArray

if TYPE_CHECKING:
    from ..worlds import BurnManWorld


class BurnManLayer(PhysicsLayer):
    """

    """

    layer_class = 'burnman'

    def __init__(self, layer_name: str, layer_index: int, world: 'BurnManWorld', burnman_layer: burnman.Layer,
                 layer_config: dict, initialize: bool = True):

        super().__init__(layer_name, layer_index, world, layer_config, initialize)

        # BurnManLayer is built on the BurnMan API and layer system.
        #    Pull out BurnMan parameters
        self._bm_layer = burnman_layer
        self._bm_material = burnman_layer.material

        # Other state variables that will be initialized later
        self._bm_mid_index = None
        self._persistent_pressure = None
        self.interp_func = None

        # Material Properties set by BurnMan
        self.interp_temperature_range = np.linspace(*tuple(self.config['interp_temperature_range']),
                                                    configurations['burnman_interpolation_N'])
        self._interp_prop_data_lookup = dict()

        # Build lookup tables
        self._build_material_property_interpolation()


    def reinit(self, initial_init: bool = False, called_from_bm_layer: bool = True):

        if initial_init:
            # Setup layer's physical properties and geometry based on the BurnMan results
            bm_radius = np.max(self.bm_layer.radii)
            bm_thickness = bm_radius - np.min(self.bm_layer.radii)
            bm_mass = self.bm_layer.mass
            self.set_geometry(radius=bm_radius, mass=bm_mass, thickness=bm_thickness)

            # Find which slice index is nearest to the middle of the layer (as determined by the radius/thickness)
            self._bm_mid_index = find_nearest(self.bm_layer.radii, self.radius - (self.thickness / 2.))

            # For BurnManLayer, assume everything is calculated by BurnMan.
            #     Override some of the items set by the set_geometry method.
            if configurations['burnman_interpolation_method'] == 'mid':
                self._pressure = self.bm_layer.pressures[self._bm_mid_index]
                self._density = self.bm_layer.density[self._bm_mid_index]
                self._gravity = self.bm_layer.gravity[self._bm_mid_index]
                self.interp_func = lambda array: array[self._bm_mid_index]
            elif configurations['burnman_interpolation_method'] == 'avg':
                self._pressure = np.average(self.bm_layer.pressures)
                self._density = np.average(self.bm_layer.density)
                self._gravity = np.average(self.bm_layer.gravity)
                self.interp_func = np.average
            elif configurations['burnman_interpolation_method'] == 'median':
                self._pressure = np.median(self.bm_layer.pressures)
                self._density = np.median(self.bm_layer.density)
                self._gravity = np.median(self.bm_layer.gravity)
                self.interp_func = np.median
            else:
                raise UnknownTidalPyConfigValue

            # Set up a pressure that will persist if the layer's state is cleared
            self._persistent_pressure = self.pressure

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

        # Call the standard reinit for the physics layer
        del called_from_bm_layer
        super().reinit(initial_init, called_from_bm_layer=True)

    def set_state(self, temperature: FloatArray = None, pressure: FloatArray = None,
                  force_update_strength: bool = True):
        """ Set the layer's temperature and update all related properties.

        Parameters
        ----------
        temperature : FloatArray = None
            Temperature of layer [K]
        pressure : FloatArray = None
            Pressure of layer [Pa]
        force_update_strength : bool = True
            If True, <layer>.update_strength() will be called

        """

        # TODO: Currently changing the temperature (or viscosity/shear) does not call BurnMan to update the full layer's equation of state.
        #    This generally should be okay for small changes, but the larger the delta-temp the larger the error will be.
        #    Recalculating the EOS can be computationally expensive, especially once it is wrapped within a bunch of TidalPy methods.
        #    This may be a major limitation to TidalPy's future development. It certainly is missing in the current release.
        #    Workaround: the user can access the burnman layer and planet to manually update the EOS, but changes won't propagate into the TidalPy classes (as of writing).

        # Call the base version first but do not let it update the strength yet.
        super().set_state(temperature, pressure, force_update_strength=False)

        if temperature is not None:
            # Set new material properties based on BurnMan Interpolation
            for prop_name, prop_data in self._interp_prop_data_lookup.items():
                setattr(self, prop_name, np.interp(self._temperature, self.interp_temperature_range, prop_data))
            self.thermal_diffusivity = self.thermal_conductivity / (self.density * self.specific_heat)

        if force_update_strength:
            # Temperature and pressure will change the strength of the layer and all of its dependencies
            self.update_strength()

    def geotherm(self, avg_temperature: float = None):
        """ Calculates layer's geotherm based on an average temperature (or layer's current temperature)

        Returns
        -------
        temperature_profile : np.ndarray
            numpy array of the adiabatic temperature profile
        """

        if avg_temperature is None:
            avg_temperature = self.temperature

        temperature_profile = burnman.geotherm.adiabatic(self.pressure_slices, avg_temperature, self.bm_material)

        return temperature_profile

    def _build_material_property_interpolation(self):
        """ Interpolates material properties based on a fixed pressure and a suggested temperature range.

        Using the BurnMan package's equation of states, this function will build interpolated lookup tables for several
            material properties. These lookup tables are functions of only temperature; pressure is assumed to not
            change significantly in a TidalPy simulation. However, the machinery to change properties based on
            pressure is built into BurnMan and can be accessed via the layer's reference to the Burnman layer:
                self.bm_material.evaluate

        The specific constant pressure used in building the lookup table is set by the layer configuration. See the
            variable 'burnman_interpolation_method' in TidalPy's main configurations.py

        Conversions between BurnMan and TidalPy property names are also performed here. Some conversions require a
             transition from molar to specific.

        The final interpolated lookup table is stored in the self._interp_prop_data_lookup which is used in the
             temperature.setter
        """

        if self.pressure is None:
            raise AttributeNotSetError

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
    def deriv_temperature(self) -> np.ndarray:
        return self._deriv_temperature

    @deriv_temperature.setter
    def deriv_temperature(self, value):
        raise ImproperPropertyHandling


    # Inner-scope properties
    # # BurnMan Layer Class
    @property
    def pressure_outer(self) -> float:
        return self._bm_layer.pressures[-1]

    @pressure_outer.setter
    def pressure_outer(self, value):
        raise ImproperPropertyHandling

    @property
    def pressure_inner(self) -> float:
        return self._bm_layer.pressures[0]

    @pressure_inner.setter
    def pressure_inner(self, value):
        raise ImproperPropertyHandling

    @property
    def temperature_outer(self) -> float:
        return self._bm_layer.temperatures[-1]

    @temperature_outer.setter
    def temperature_outer(self, value):
        raise ImproperPropertyHandling

    @property
    def temperature_inner(self) -> float:
        return self._bm_layer.temperatures[0]

    @temperature_inner.setter
    def temperature_inner(self, value):
        raise ImproperPropertyHandling

    @property
    def density_outer(self) -> float:
        return self._bm_layer.density[-1]

    @density_outer.setter
    def density_outer(self, value):
        raise ImproperPropertyHandling

    @property
    def density(self) -> float:
        return self._density

    @density.setter
    def density(self, value: float):
        if type(value) is not float:
            raise IncorrectAttributeType
        self._density = value

    @property
    def density_inner(self) -> float:
        return self._bm_layer.density[0]

    @density_inner.setter
    def density_inner(self, value):
        raise ImproperPropertyHandling

    @property
    def gravity_outer(self) -> float:
        return self._bm_layer.gravity[-1]

    @gravity_outer.setter
    def gravity_outer(self, value):
        raise ImproperPropertyHandling

    @property
    def gravity(self) -> float:
        return self._gravity

    @gravity.setter
    def gravity(self, value: float):
        if type(value) is not float:
            raise IncorrectAttributeType
        self._gravity = value

    @property
    def gravity_inner(self) -> float:
        return self._bm_layer.gravity[0]

    @gravity_inner.setter
    def gravity_inner(self, value):
        raise ImproperPropertyHandling

    @property
    def radii(self) -> np.ndarray:
        return self._bm_layer.radii

    @radii.setter
    def radii(self, value):
        raise ImproperPropertyHandling

    # TODO: make a "slice" inner class that handles this stuff an 3d Love calculation?
    @property
    def pressure_slices(self) -> np.ndarray:
        return self._bm_layer.pressures

    @pressure_slices.setter
    def pressure_slices(self, value):
        raise ImproperPropertyHandling

    @property
    def density_slices(self) -> np.ndarray:
        return self._bm_layer.density

    @density_slices.setter
    def density_slices(self, value):
        raise ImproperPropertyHandling

    @property
    def gravity_slices(self) -> np.ndarray:
        return self._bm_layer.gravity

    @gravity_slices.setter
    def gravity_slices(self, value):
        raise ImproperPropertyHandling


    # Aliased properties
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
