from __future__ import annotations

import copy
from typing import TYPE_CHECKING, Tuple, Union, Dict, List

import burnman
import numpy as np

from ...utilities.arrayHelp.shaper import reshape_help
from ...tides.tides import Tides
from .defaults import layer_defaults
from ... import debug_mode
from ...burnman_interface.conversion import burnman_property_name_conversion, burnman_property_value_conversion
from ...configurations import burnman_interpolation_N, burnman_interpolation_method
from ...exceptions import (AttributeNotSetError, ImproperAttributeHandling, IncorrectAttributeType,
                           ParameterMissingError, ReinitError, UnknownTidalPyConfigValue, UnusualRealValueError,
                           OuterscopeAttributeSetError)
from ...initialize import log
from ...radiogenics.radiogenics import Radiogenics
from ...structures.physical import PhysicalObjSpherical
from ...cooling import CoolingModel, CoolingOutputTypeArray
from ...rheology import Rheology
from ...types import floatarray_like, NoneType
from ...utilities.dictionary_utils import nested_get
from ...utilities.numpy_help import find_nearest
from ...types import ArrayNone, FloatArray


if TYPE_CHECKING:
    from ..worlds import ThermalWorld


class ThermalLayer(PhysicalObjSpherical):
    """ Layer Class: Tidal Planets are made up of various Layers of different compositions.

    The ThermalLayer class acts as the holder for:
        - Material properties (initially set by configurations or BurnMan data)
        - Rheology Class
            * ComplexCompliance Class
            * Viscosity Class
            * PartialMelting Class
        - Cooling Class
        - Radiogenic Class

    """

    default_config = layer_defaults
    layer_class = 'thermal'

    def __init__(self, layer_name: str, layer_index: int, world: 'ThermalWorld', burnman_layer: burnman.Layer,
                 layer_config: dict, initialize: bool = True):

        # Load layer defaults based on layer type
        self.type = layer_config['type']
        self.default_config = self.default_config[self.type]

        # Setup Physical Layer Geometry
        super().__init__(layer_config)

        # State properties
        self._layer_index = layer_index
        self._world = world
        # Pull out information from the already initialized burnman Layer
        self._bm_layer = burnman_layer
        self._bm_material = burnman_layer.material
        # Setup by burnman interpolations
        self._gravity = None
        self._density = None
        # Set later on by user
        self._pressure = None
        self._temperature = None
        # State Derivatives
        self._deriv_temperature = None  # type: ArrayNone

        # Other attributes
        self.name = layer_name
        self.tidal_scale = 1.
        self.heat_sources = None
        # Flags
        self.is_top_layer = layer_index == 0
        self.is_tidal = False
        self.use_pressure_in_strength_calc = True

        # Model holders - setup in reinit()
        self.rheology = None       # type: Union['Rheology', NoneType]
        self.cooling_model = None  # type: Union['CoolingModel', NoneType]
        self.radiogenics = None    # type: Union['Radiogenics', NoneType]

        # Setup geometry based on the BurnMan results
        bm_radius = np.max(self.bm_layer.radii)
        bm_thickness = bm_radius - np.min(self.bm_layer.radii)
        bm_mass = self.bm_layer.mass
        self.set_geometry(radius=bm_radius, mass=bm_mass, thickness=bm_thickness)
        self._bm_mid_index = find_nearest(self.bm_layer.radii, self.radius - self.thickness / 2.)

        # Attributes calculated by BurnMan but set at a specific spot
        if burnman_interpolation_method == 'mid':
            self.pressure = self.bm_layer.pressures[self._bm_mid_index]
            self.density = self.bm_layer.density[self._bm_mid_index]
            self.gravity = self.bm_layer.gravity[self._bm_mid_index]
            self.interp_func = lambda array: array[self._bm_mid_index]
        elif burnman_interpolation_method == 'avg':
            self.pressure = np.average(self.bm_layer.pressures)
            self.density = np.average(self.bm_layer.density)
            self.gravity = np.average(self.bm_layer.gravity)
            self.interp_func = np.average
        elif burnman_interpolation_method == 'median':
            self.pressure = np.median(self.bm_layer.pressures)
            self.density = np.median(self.bm_layer.density)
            self.gravity = np.median(self.bm_layer.gravity)
            self.interp_func = np.median
        else:
            raise UnknownTidalPyConfigValue

        # Setup Physical Layer Geometry
        self.slices = self.config['slices'] # TODO: ????
        self.material_name = self.config['material']



        # TODO: ????
        # Material Properties set by BurnMan
        self.bulk_modulus = None
        self.thermal_expansion = None
        self.specific_heat = None
        self.energy_per_therm = None
        self.interp_temperature_range = np.linspace(*tuple(self.config['interp_temperature_range']),
                                                    burnman_interpolation_N)
        self._interp_prop_data_lookup = dict()
        self._build_material_property_interpolation()

        # Material properties that might have been affected by new configuration file (BurnMan does not calculate these)
        self.static_shear_modulus = None  # type: ArrayNone
        self.thermal_diffusivity = None  # type: ArrayNone
        self.thermal_conductivity = None  # type: ArrayNone
        self.heat_fusion = None  # type: ArrayNone
        self.temp_ratio = None  # type: ArrayNone
        self.stefan = None  # type: ArrayNone

        if initialize:
            self.reinit()

    def reinit(self):

        super().reinit()

        # Load in configurations
        if self.config['use_tvf']:
            self.tidal_scale = self.volume / self.world.volume
        self.is_tidal = self.config['is_tidally_active']
        self.use_pressure_in_strength_calc = self.config['use_pressure_in_strength_calc']

        if self.config['use_surf_gravity']:
            # Use surface gravity for layer instead of the gravity set by interpolating burnman data (mid/avg/etc)
            # This primarily affects convection calculation
            self.gravity = self.gravity_surf

        # Material properties that might have been affected by new configuration files
        self.static_shear_modulus = self.config['shear_modulus']
        self.heat_fusion = self.config['heat_fusion']
        self.thermal_conductivity = self.config['thermal_conductivity']
        self.temp_ratio = self.config['boundary_temperature_ratio']
        self.stefan = self.config['stefan']

        # Setup Models
        self.rheology = Rheology(self, store_config_in_layer=True)
        self.radiogenics = Radiogenics(self, store_config_in_layer=True)
        self.cooling_model = CoolingModel(self, store_config_in_layer=True)

        # Find all heating sources
        self.heat_sources = list()
        if 'off' not in self.radiogenics.model:
            self.heat_sources.append(lambda : self.radiogenic_heating)
        if self.is_tidal:
            self.heat_sources.append(lambda: self.tidal_heating)

    def set_state(self, temperature: FloatArray = None, pressure: FloatArray = None):
        """ Set the layer's temperature and update all related properties.

        Parameters
        ----------
        temperature : FloatArray = None
            Temperature of layer [K]
        pressure : FloatArray = None
            Pressure of layer [Pa]

        """

        if temperature is not None:
            temperature = reshape_help(temperature, self.world.global_shape)
            self._temperature = temperature
        if pressure is not None:
            pressure = reshape_help(pressure, self.world.global_shape)
            self._pressure = pressure

        # Temperature and pressure will change the strength of the layer and all of its dependencies
        self.update_strength()

    def set_temperature(self, temperature: FloatArray):
        """ Wrapper for ThermalLayer.set_state

        Parameters
        ----------
        temperature : FloatArray
            Temperature of layer [K]

        See Also
        --------
        TidalPy.structures.layers.layers.ThermalLayer.set_state
        """

        self.set_state(temperature=temperature, pressure=None)

    def set_pressure(self, pressure: FloatArray):
        """ Wrapper for ThermalLayer.set_state

        Parameters
        ----------
        pressure : FloatArray
            Pressure of layer [Pa]

        See Also
        --------
        TidalPy.structures.layers.layers.ThermalLayer.set_state
        """

        self.set_state(temperature=None, pressure=pressure)

    def set_strength(self, viscosity: FloatArray = None, shear_modulus: FloatArray = None):
        """ Manual set the viscosity and shear modulus of the layer, independent of temperature.

        This method by-passes the self.viscosity_func and allows the user to manually set the viscosity and shear
        of the layer.

        Future : TODO
        ------
        Currently it does not change the layer's temperature. In the future an effective temperature can be
            (optionally) calculated. Make sure to implement this at the Rheology class level.

        Parameters
        ----------
        viscosity : FloatArray
            The new viscosity of the layer in [Pa s]
        shear_modulus : FloatArray
            The new shear modulus of the layer in [Pa]
        """

        # Pass the new strength to the rheology class
        self.rheology.set_state(viscosity, shear_modulus, called_by_layer=True)

        self.update_strength(already_called_rheology=True)

    def update_radiogenics(self, force_calculation: bool = False) -> np.ndarray:
        """ Calculate the radiogenic heating rate within this layer, requires the world.time to be set.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.

        Returns
        -------
        radiogenic_heating : np.ndarray
            Radiogenic heating rate in [Watts]
        """

        try:
            self.radiogenics.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error
            else:
                pass

        # New radiogenic heating will change the temperature derivative
        self.calc_temperature_derivative(force_calculation=False)

        return self.radiogenic_heating

    def update_cooling(self, force_calculation: bool = False, call_deriv_calc: bool = True):
        """ Calculate parameters related to cooling of the layer, be it convection or conduction

        Using the self.cooling class and the current temperature and viscosity, this method will calculate convection
            parameters. Some of these may be zeros if convection is forced off in the planet's configuration.

        By default this method will ignore parameter missing errors to avoid initialization errors. This behaviour
            can be altered by the force_calculation flag.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.
        call_deriv_calc : bool = True
            If True, then this method will attempt to call self.calc_temperature_derivative
        """

        try:
            self.cooling_model.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error

        if self.is_top_layer and self.world.orbit is not None:
            # TODO: Add some sort of convergence here?
            # Tell the planet that there is a new surface flux and thus a new surface temperature
            self.world.update_surface_temperature()

        if call_deriv_calc:
            # New cooling will change the temperature derivative
            self.calc_temperature_derivative(force_calculation=False)

    def update_tides(self, force_calculation: bool = False, call_deriv_calc: bool = True):
        """ Calculate, and update the layer's state variables with tidal heating and torques based on its current state

        Using the self.tides class and the layer's current temperature, viscosity, shear modulus, and frequencies;
            the layer will attempt to calculate a new effective rigidity, complex compliance, and complex love number,
            ultimately it will find tidal heating and torque then inform the planet that it should update its own
            global calculations.

        Viscosity, Shear Modulus, and Orbital Freq are required inputs to most rheological models. If the planet is not
            forced into spin-sync, then Spin Freq is also required. If a missing parameter is raised then the method
            will exit unless force_calculation is set to True.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.
        call_deriv_calc : bool = True
            If True, then this method will attempt to call self.calc_temperature_derivative

        """

        if self.is_tidal:
            try:
                self.world.update_tides()
            except ParameterMissingError as error:
                if force_calculation:
                    raise error

        if call_deriv_calc:
            self.calc_temperature_derivative(force_calculation=False)

    def update_strength(self, already_called_rheology: bool = False):
        """ Calculate viscosity and shear modulus and update the layer's state with new values.

        Rheology class handles the change in pre-melt shear modulus and viscosity as well as the impact of partial
            melting. It will then calculate the complex compliances which will impact the global tides.

        The pressure of the layer will be used, if present, in the calculations (assuming that the viscosity and shear
            functions actually depend upon pressure). Otherwise, a pressure of zero will be assumed.

        Parameters
        ----------
        already_called_rheology : bool = False
            Some other methods of this class may call this method and may have already performed some of the tasks
                that this method would normally perform.

        """

        if not already_called_rheology:
            # Tell rheology class to update strength (calculation of effective viscosity, shear, and complex compliance)
            self.rheology.update_strength(called_from_layer=True)

        # Changing the strength will change the cooling properties. Don't have it call derivative calculation yet.
        self.update_cooling(call_deriv_calc=False)

        # Changing the strength will change tides
        self.update_tides(call_deriv_calc=False)

        # The above changes will impact the temperature derivative calculation. So let's call that now.
        self.calc_temperature_derivative(force_calculation=False)

    def calc_temperature_derivative(self, force_calculation: bool = False) -> np.ndarray:
        """ Calculate the change in temperature within this layer over time.

        dT/dt = (heating_in - cooling_out) / (M * c_p * (St + 1) * Tr)

        Returns
        -------
        dT/dt : np.ndarray
            The derivative of temperature with respect to time.
        """

        # Initialize total heating based on if there is a layer warming this one from below
        if self.layer_below is None:
            total_heating = np.zeros_like(self.temperature)
        elif self.layer_below.heat_flux is None:
            total_heating = np.zeros_like(self.temperature)
        else:
            total_heating = self.layer_below.heat_flux * self.surface_area_inner

        try:
            # Heat sources are added in the reinit() method.
            for heat_source_func in self.heat_sources:
                heat_source = heat_source_func()
                if heat_source is None:
                    raise ParameterMissingError
                total_heating += heat_source

            # Cooling is calculated by the cooling_model class
            total_cooling = self.cooling

            self._deriv_temperature = (total_heating - total_cooling) / \
                                     (self.mass * self.specific_heat * (self.stefan + 1.) * self.temp_ratio)
        except ParameterMissingError as error:
            if force_calculation:
                raise error

        return self.deriv_temperature

    def geotherm(self, avg_temperature: float = None):
        """ Calculates layer's geotherm

        Returns
        -------
        temperature_profile : np.ndarray
            numpy array of the adiabatic temperature profile
        """

        if avg_temperature is None:
            avg_temperature = self.temperature

        temperature_profile = burnman.geotherm.adiabatic(self.pressure_slices, avg_temperature, self.bm_material)

        return temperature_profile

    # State properties
    @property
    def layer_index(self) -> int:
        return self._layer_index

    @layer_index.setter
    def layer_index(self, value):
        raise ImproperAttributeHandling('Layer index can not be changed after planet has been initialized.')

    @property
    def world(self) -> 'ThermalWorld':
        return self._world

    @world.setter
    def world(self, value):
        raise ImproperAttributeHandling('Can not change world association after a layer has been initialized.')

    @property
    def bm_layer(self) -> burnman.Layer:
        return self._bm_layer

    @bm_layer.setter
    def bm_layer(self, value):
        raise ImproperAttributeHandling('Can not change burnman layer information after layer has been initialized.')

    @property
    def bm_material(self) -> burnman.Material:
        return self._bm_material

    @bm_material.setter
    def bm_material(self, value):
        raise ImproperAttributeHandling('Can not change burnman layer information after layer has been initialized.')

    @property
    def pressure_upper(self) -> float:
        return self._bm_layer.pressures[-1]

    @pressure_upper.setter
    def pressure_upper(self, value):
        raise ImproperAttributeHandling

    @property
    def pressure(self) -> np.ndarray:
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        self.set_pressure(value)

    @property
    def pressure_lower(self) -> float:
        return self._bm_layer.pressures[0]

    @pressure_lower.setter
    def pressure_lower(self, value):
        raise ImproperAttributeHandling

    @property
    def temperature_upper(self) -> float:
        return self._bm_layer.temperatures[-1]

    @temperature_upper.setter
    def temperature_upper(self, value):
        raise ImproperAttributeHandling

    @property
    def temperature(self) -> np.ndarray:
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self.set_temperature(value)

    @property
    def temperature_lower(self) -> float:
        return self._bm_layer.temperatures[0]

    @temperature_lower.setter
    def temperature_lower(self, value):
        raise ImproperAttributeHandling

    @property
    def density_upper(self) -> float:
        return self._bm_layer.density[-1]

    @density_upper.setter
    def density_upper(self, value):
        raise ImproperAttributeHandling

    @property
    def density(self) -> float:
        return self._density

    @density.setter
    def density(self, value: float):
        if type(value) is not float:
            raise IncorrectAttributeType
        self._density = value

    @property
    def density_lower(self) -> float:
        return self._bm_layer.density[0]

    @density_lower.setter
    def density_lower(self, value):
        raise ImproperAttributeHandling

    @property
    def gravity_upper(self) -> float:
        return self._bm_layer.gravity[-1]

    @gravity_upper.setter
    def gravity_upper(self, value):
        raise ImproperAttributeHandling

    @property
    def gravity(self) -> float:
        return self._gravity

    @gravity.setter
    def gravity(self, value: float):
        if type(value) is not float:
            raise IncorrectAttributeType
        self._gravity = value

    @property
    def gravity_lower(self) -> float:
        return self._bm_layer.gravity[0]

    @gravity_lower.setter
    def gravity_lower(self, value):
        raise ImproperAttributeHandling

    @property
    def radii(self) -> np.ndarray:
        return self._bm_layer.radii

    @radii.setter
    def radii(self, value):
        raise ImproperAttributeHandling

    # TODO: make a "slice" inner class that handles this stuff an 3d Love calculation?
    @property
    def pressure_slices(self) -> np.ndarray:
        return self._bm_layer.pressures

    @pressure_slices.setter
    def pressure_slices(self, value):
        raise ImproperAttributeHandling

    @property
    def density_slices(self) -> np.ndarray:
        return self._bm_layer.density

    @density_slices.setter
    def density_slices(self, value):
        raise ImproperAttributeHandling

    @property
    def gravity_slices(self) -> np.ndarray:
        return self._bm_layer.gravity

    @gravity_slices.setter
    def gravity_slices(self, value):
        raise ImproperAttributeHandling

    @property
    def deriv_temperature(self) -> np.ndarray:
        return self._deriv_temperature

    @deriv_temperature.setter
    def deriv_temperature(self, value):
        raise ImproperAttributeHandling


    # Inner-scope properties
    # # Rheology Class
    @property
    def viscosity(self):
        return self.rheology.postmelt_viscosity

    @viscosity.setter
    def viscosity(self, value):
        self.rheology.postmelt_viscosity = value

    @property
    def liquid_viscosity(self):
        return self.rheology.liquid_viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        self.rheology.liquid_viscosity = value

    @property
    def shear_modulus(self):
        return self.rheology.postmelt_shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        self.rheology.postmelt_shear_modulus = value

    @property
    def compliance(self):
        return self.rheology.postmelt_compliance

    @compliance.setter
    def compliance(self, value):
        self.rheology.postmelt_compliance = value

    @property
    def complex_compliances(self):
        return self.rheology.complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        self.rheology.complex_compliances = value

    @property
    def melt_fraction(self):
        return self.rheology.melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        self.rheology.melt_fraction = value

    # # Radiogenics Class
    @property
    def radiogenic_heating(self):
        return self.radiogenics.heating

    @radiogenic_heating.setter
    def radiogenic_heating(self, value):
        self.radiogenics.heating = value

    # # Cooling Class
    @property
    def cooling(self):
        return self.cooling_model.cooling

    @cooling.setter
    def cooling(self, value):
        self.cooling_model.cooling = value

    @property
    def cooling_flux(self):
        return self.cooling_model.cooling_flux

    @cooling_flux.setter
    def cooling_flux(self, value):
        self.cooling_model.cooling_flux = value

    @property
    def boundary_layer_thickness(self):
        return self.cooling_model.boundary_layer_thickness

    @boundary_layer_thickness.setter
    def boundary_layer_thickness(self, value):
        self.cooling_model.boundary_layer_thickness = value

    @property
    def rayleigh(self) -> np.ndarray:
        return self.cooling_model.rayleigh

    @rayleigh.setter
    def rayleigh(self, value):
        self.cooling_model.rayleigh = value

    @property
    def nusselt(self) -> np.ndarray:
        return self.cooling_model.nusselt

    @nusselt.setter
    def nusselt(self, value):
        self.cooling_model.nusselt = value


    # Outer-scope properties
    # # World Class
    @property
    def layer_below(self) -> Union['ThermalLayer', NoneType]:
        if self.layer_index == 0:
            return None
        else:
            return self.world.layers_by_index[self.layer_index - 1]

    @layer_below.setter
    def layer_below(self, value):
        raise OuterscopeAttributeSetError

    @property
    def layer_above(self) -> Union['ThermalLayer', NoneType]:
        if self.layer_index == self.world.num_layers - 1:
            return None
        else:
            return self.world.layers_by_index[self.layer_index + 1]

    @layer_above.setter
    def layer_above(self, value):
        raise OuterscopeAttributeSetError

    @property
    def time(self):
        return self.world.time

    @time.setter
    def time(self, value):
        raise OuterscopeAttributeSetError

    # # Tides Class
    @property
    def tidal_heating(self):
        return self.world.tides.tidal_heating_by_layer[self]

    @tidal_heating.setter
    def tidal_heating(self, value):
        raise OuterscopeAttributeSetError


    # Aliased properties
    @property
    def gravity_surf(self):
        return self.gravity_upper

    @gravity_surf.setter
    def gravity_surf(self, value):
        self.gravity_upper = value

    @property
    def blt(self):
        return self.boundary_layer_thickness

    @blt.setter
    def blt(self, value):
        self.boundary_layer_thickness = value


    # Dunder methods
    def __str__(self):

        if self.world is None:
            text = f'[Layer {self.name.title()}:{self.type.title()} no world]'
        else:
            text = f'[Layer {self.name.title()}:{self.type.title()} in {self.world}]'
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















# TODO OLD

    def _build_material_property_interpolation(self):
        """ Interpolates material properties based on a fixed pressure and a suggested temperature range.

        Using the Burnman package's equation of states, this function will build interpolated lookup tables for several
        material properties. These lookup tables are functions of only temperature; pressure is assumed to not change
        significantly in a TidalPy simulation. However, the machinery to change properties based on temperature is
        built into Burnman and can be accessed via the layer's reference to the Burnman layer:
        self.bm_material.evaluate

        The specific constant pressure used in building the lookup table is set by the layer configuration. See the
        variable 'burnman_interpolation_method' in the layers __init__ method.

        Conversions between Burnman TidalPy property names are also performed here. Some conversions require a
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



    def set_geometry(self, radius: float, mass: float, thickness: float = None):

        super().set_geometry(radius, mass, thickness)




    @property
    def temperature(self) -> np.ndarray:

        return self._temperature

    @temperature.setter
    def temperature(self, value):

        if debug_mode:
            if type(value) not in floatarray_like:
                raise IncorrectAttributeType
            if type(value) == np.ndarray:
                if np.any(value > 1.0e5) or np.any(value < 5.0):
                    raise UnusualRealValueError
            else:
                if value > 1.0e5 or value < 5.0:
                    raise UnusualRealValueError

        self._temperature = value_cleanup(value)

        # Update material properties and phase changes
        self._melt_fraction = self.partial_melt.calc_melt_fraction()
        # Set new material properties based on BurnMan Interpolation
        for prop_name, prop_data in self._interp_prop_data_lookup.items():
            setattr(self, prop_name, np.interp(self._temperature, self.interp_temperature_range, prop_data))
        self.thermal_diffusivity = self.thermal_conductivity / (self.density * self.specific_heat)
        # Update shear and viscosity
        self.update_strength()

        # Update other models
        # TODO: Are these forces still needed or are the a relic of past code?
        force_cooling = True
        force_tides = False
        if self.orbital_freq is not None:
            if self.spin_freq is not None or self.is_spin_sync:
                force_tides = True

        self.update_tides(force_calculation=force_tides)
        self.update_cooling(force_calculation=force_cooling)

    # TODO: For now lower temperature just returns self.temperature. This is used in cooling calculations
    @property
    def temperature_lower(self):
        return self.temperature

    @temperature_lower.setter
    def temperature_lower(self, value):
        raise ImproperAttributeHandling('temperature_lower is calculated by setting self.temperature.')

    @property
    def temperature_surf(self):

        if self.is_top_layer:
            if self.world.orbit is None:
                # surface temperature will not have been set yet. Use a fake value for now
                log('Surface temperature is not set until an orbit is applied to the planet. '
                    'Using default value of 100K.', level='info')
                return np.asarray([150.])
            else:
                return self.world.surface_temperature

        if self.layer_above is not None:
            if self.layer_above.temperature_lower is not None:
                return self.layer_above.temperature_lower

        raise ParameterMissingError
