from __future__ import annotations

import copy
from typing import TYPE_CHECKING, Tuple

import burnman
import numpy as np

from TidalPy.tides.tides import Tides
from .defaults import layer_defaults
from .. import debug_mode
from ..bm.conversion import burnman_property_name_conversion, burnman_property_value_conversion
from ..configurations import burnman_interpolation_N, burnman_interpolation_method
from ..exceptions import (AttributeNotSetError, ImproperAttributeHandling, IncorrectAttributeType,
                          ParameterMissingError, ReinitError, UnknownTidalPyConfigValue, UnusualRealValueError)
from ..initialize import log
from ..radiogenics.radiogenics import Radiogenics
from ..structures.physical import PhysicalObjSpherical
from ..thermal import find_viscosity
from ..thermal.cooling import Cooling
from ..thermal.partial_melt import PartialMelt
from ..types import floatarray_like
from ..utilities.dict_tools import nested_get
from ..utilities.numpy_help import find_nearest, value_cleanup


if TYPE_CHECKING:
    from .worlds import TidalWorld


class ThermalLayer(PhysicalObjSpherical):
    """ Layer Class: Tidal Planets are made up of various Layers of different compositions.

        Any functionality that requires information about material properties should be called from the layer class.
    """

    default_config = copy.deepcopy(layer_defaults)
    class_type = 'thermal'

    def __init__(self, layer_name: str, world: TidalWorld, burnman_layer: burnman.Layer, layer_config: dict):

        # Load layer defaults
        self.type = layer_config['type']
        self.default_config = self.default_config[self.type]

        # Setup Physical Layer Geometry
        super().__init__(layer_config, call_reinit=False)

        self.name = layer_name
        self.world = world

        # Information about the other layers (mostly for gravity and cooling calculations)
        self.layer_below = None  # type: ThermalLayer or None
        self.layer_above = None  # type: ThermalLayer or None

        # Pull out information from the already initialized burnman Layer
        self.bm_layer = burnman_layer
        self.bm_material = burnman_layer.material

        # Other Attributes - these should not change on a reinit
        self.pressure_upper = self.bm_layer.pressures[-1]
        self.pressure_lower = self.bm_layer.pressures[0]
        self.density_upper = self.bm_layer.density[-1]
        self.density_lower = self.bm_layer.density[0]

        # Setup geometry based on the BurnMan results
        bm_radius = np.max(self.bm_layer.radii)
        bm_thickness = bm_radius - np.min(self.bm_layer.radii)
        bm_mass = self.bm_layer.mass
        self.set_geometry(radius=bm_radius, mass=bm_mass, thickness=bm_thickness)
        self.bm_mid_index = find_nearest(self.bm_layer.radii, self.radius - self.thickness / 2.)
        self._gravity_surf = self.bm_layer.gravity[-1]

        # Pull out slice information
        self.radii = self.bm_layer.radii
        self.pressure_slices = self.bm_layer.pressures
        self.density_slices = self.bm_layer.density
        self.gravity_slices = self.bm_layer.gravity

        # Attributes calculated by BurnMan but set at a specific spot
        if burnman_interpolation_method == 'mid':
            self.pressure = self.bm_layer.pressures[self.bm_mid_index]
            self.density = self.bm_layer.density[self.bm_mid_index]
            self.gravity = self.bm_layer.gravity[self.bm_mid_index]
            self.interp_func = lambda array: array[self.bm_mid_index]
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
        self.slices = self.config['slices']
        self.material_name = self.config['material']

        # Load in switches
        if self.config['use_surf_gravity']:
            # This primarily affects convection calculation
            self.gravity = self.gravity_surf

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
        self.static_shear_modulus = None  # type: float or None
        self.thermal_diffusivity = None  # type: float or None
        self.thermal_conductivity = None  # type: float or None
        self.heat_fusion = None  # type: float or None
        self.temp_ratio = None  # type: float or None
        self.stefan = None  # type: float or None

        # State Variables
        self._temperature = None  # type: np.ndarray or None
        self._temperature_lower = None  # type: np.ndarray or None
        self._melt_fraction = None  # type: np.ndarray or None
        self._viscosity = None  # type: np.ndarray or None
        self._shear_modulus = None  # type: np.ndarray or None
        self._heat_flux = None  # type: np.ndarray or None
        self._rayleigh = None  # type: np.ndarray or None
        self._nusselt = None  # type: np.ndarray or None
        self._blt = None  # type: np.ndarray or None
        self._effective_rigidity = None  # type: np.ndarray or None
        self._complex_compliance = None  # type: np.ndarray or None
        self._complex_love = None  # type: np.ndarray or None
        self._tidal_heating = None  # type: np.ndarray or None
        self._radiogenic_heating = None  # type: np.ndarray or None

        # State Derivatives
        self._diff_temperature = None  # type: np.ndarray or None

        # Model Holders
        self.cooling = None  # type: Cooling or None
        self.partial_melt = None  # type: PartialMelt or None
        self.radiogenics = None  # type: Radiogenics or None
        self.rheology = None  # type: Tides or None
        self.heat_sources = None

        # Function Holders
        self.viscosity_func, self.viscosity_inputs, self.viscosity_live_inputs = None, None, None
        self.viscosity_liq_func, self.viscosity_liq_inputs, self.viscosity_liq_live_inputs = None, None, None

        # Constants
        self.tidal_scale = 1.
        if self.config['use_tvf']:
            self.tidal_scale = self.volume / self.world.volume

        # Flags
        self.is_top_layer = False
        self.is_tidal = False

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

    def reinit(self):

        super().reinit()
        log(f'Reinit called for Layer: {self.name} in {self.world.name}.', level='debug')

        # Check for changes to config that may break the planet
        if self._old_config is not None:
            for critical_attribute in ['material', 'material_source', 'type', 'thickness', 'radius', 'slices']:
                if self.config[critical_attribute] != self._old_config[critical_attribute]:
                    raise ReinitError

        # Material properties that might have been affected by new configuration files
        self.static_shear_modulus = np.asarray([self.config['shear_modulus']])
        self.heat_fusion = self.config['heat_fusion']
        self.thermal_conductivity = self.config['thermal_conductivity']
        self.temp_ratio = self.config['boundary_temperature_ratio']
        self.stefan = self.config['stefan']

        # Setup viscosity functions
        solid_visco_model = nested_get(self.config, ['rheology', 'solid_viscosity', 'model'], default=None)
        liquid_visco_model = nested_get(self.config, ['rheology', 'liquid_viscosity', 'model'], default=None)
        solid_visco_params = None
        if 'rheology' in self.config:
            if 'solid_viscosity' in self.config['rheology']:
                solid_visco_params = self.config['rheology']['solid_viscosity']
        liquid_visco_params = None
        if 'rheology' in self.config:
            if 'liquid_viscosity' in self.config['rheology']:
                liquid_visco_params = self.config['rheology']['liquid_viscosity']

        self.viscosity_func, self.viscosity_inputs, self.viscosity_live_inputs = \
            find_viscosity(solid_visco_model, parameters=solid_visco_params,
                           default_key=[self.type, 'solid_viscosity'])
        self.viscosity_liq_func, self.viscosity_liq_inputs, self.viscosity_liq_live_inputs = \
            find_viscosity(liquid_visco_model, parameters=liquid_visco_params,
                           default_key=[self.type, 'liquid_viscosity'])

        # Heat sources used in temperature derivative calculation
        self.heat_sources = list()

        # Setup Partial Melting
        self.partial_melt = PartialMelt(self)

        # Setup Cooling
        self.cooling = Cooling(self)

        # Setup Radiogenics
        self.radiogenics = Radiogenics(self)
        self.heat_sources.append(lambda: self.radiogenic_heating)

        # Setup Tides
        self.rheology = Tides(self)
        if self.rheology.model != 'off':
            self.is_tidal = True
            self.heat_sources.append(lambda: self.tidal_heating)
        else:
            self.is_tidal = False

        self.set_geometry(self.radius, self.mass, self.thickness)

    def set_geometry(self, radius: float, mass: float, thickness: float = None):

        super().set_geometry(radius, mass, thickness)

    def set_strength(self, viscosity: np.ndarray = None, shear_modulus: np.ndarray = None):
        """ Manual set the viscosity and shear modulus of the layer, independent of temperature.

        This method by-passes the self.viscosity_func and allows the user to manually set the viscosity and shear
        of the layer.

        This method reproduces the functionality of the viscosity.setter and shear_modulus.setter, but allows them
        to be changed simultaneously leading to slightly better performance. Use the default setters if you are only
        changing one of the parameters at a time.

        Future
        ------
            Currently it does not change the layer's temperature. In the future an effective temperature can be
            (optionally) calculated.

        Parameters
        ----------
        viscosity : np.ndarray
            The new viscosity of the layer in [Pa s]
        shear_modulus : np.ndarray
            The new shear modulus of the layer in [Pa]
        """

        if debug_mode:
            for value in [viscosity, shear_modulus]:
                if value is None:
                    continue
                if type(value) not in floatarray_like:
                    raise IncorrectAttributeType
                if type(value) == np.ndarray:
                    if np.any(value > 1.0e30) or np.any(value < 1.0e-10):
                        raise UnusualRealValueError
                else:
                    if value > 1.0e30 or value < 1.0e-10:
                        raise UnusualRealValueError

        if viscosity is not None:
            self._viscosity = value_cleanup(viscosity)

        if shear_modulus is not None:
            self._shear_modulus = value_cleanup(shear_modulus)

        # TODO: Have an option to calculate an effective temperature given the viscosity and shear modulus.
        #   If this is implemented then make sure it makes its way to the viscosity and shear .setters

        self.update_cooling()
        self.update_tides()
        if self.world is not None and self.world.tidal_heating_coeffs is not None and \
                self.world.tidal_ztorque_coeffs is not None:
            self.world.update_global_tides()

    def update_radiogenics(self, force_calculation: bool = False):
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
            self._radiogenic_heating = self.radiogenics.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error
            else:
                pass

        return self.radiogenic_heating

    def update_cooling(self, force_calculation: bool = False):
        """ Calculate parameters related to cooling of the layer, be it convection or conduction

        Using the self.cooling class and the current temperature and viscosity, this method will calculate convection
        parameters. Some of these may be zeros if convection is forced off in the planet's configuration.

        By default this method will ignore parameter missing errors to avoid initialization errors. This behaviour
        can be altered by the force_calculation flag.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.

        Returns
        -------
        cooling_flux : np.ndarray
            The flux leaving this layer and entering the layer above it in [Watts m-2]
        blt : np.ndarray
            Boundary layer thickness of the layer in [m]. The max thickness is 50% of layer thickness
        rayleigh : np.ndarray
            Rayleigh number of the layer. It is equal to zero if convection is forced off or is simply very weak
        nusselt : np.ndarray
            Nusselt number of the layer. It is equal to one if convection is forced off or is very weak

        """

        try:
            cooling_flux, blt, rayleigh, nusselt = self.cooling.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error
            else:
                cooling_flux, blt, rayleigh, nusselt = None, None, None, None

        self._heat_flux, self._blt, self._rayleigh, self._nusselt = cooling_flux, blt, rayleigh, nusselt

        if self.is_top_layer and self.world.orbit is not None:
            # TODO: Add some sort of convergence here?
            # Tell the planet that there is a new surface flux and thus a new surface temperature
            self.world.update_surface_temperature()

        return cooling_flux, blt, rayleigh, nusselt

    def update_tides(self, force_calculation: bool = False, call_from_world: bool = False) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Calculate, and update the layer's state with, tidal and rheological parameters based on current state

        Using the self.rheology class and the layer's current temperature, viscosity, shear modulus, and tidal modes,' \
        the layer will attempt to calculate a new effective rigidity, complex compliance, complex love number, and
        tell the planet that it should update its own global calculations.

        Viscosity, Shear Modulus, and Tidal Modes are required inputs to most rheological models. If a missing parameter
        is raised then the method will exit unless force_calculation is set to True.

        Parameters
        ----------
        force_calculation : bool
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.

        Returns
        -------
        effective_rigid : np.ndarray or None
            Effective rigidity of the layer (dependent upon shear modulus) [unitless]
            The order can be changed via configuration files, default is 2nd order.
        complex_comp : np.ndarray or None
            Complex compliance of the layer in [Pa-1]
        complex_love : np.ndarray or None
            Complex Love number, k, of the layer in [Pa-1]
            The order can be changed via configuration files, default is 2nd order.
        """

        if self.tidal_modes is None:
            num_modes = 1
        else:
            num_modes = len(self.tidal_modes)

        update_world_tides = self.world.orbit is not None
        try:
            effective_rigid, complex_comp, complex_love = self.rheology.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error
            else:
                update_world_tides = False
                effective_rigid = tuple([np.asarray(0., dtype=np.float) for _ in range(num_modes)])
                complex_comp = tuple([np.asarray(0., dtype=np.complex) for _ in range(num_modes)])
                complex_love = tuple([np.asarray(0., dtype=np.complex) for _ in range(num_modes)])

        self._effective_rigidity, self._complex_compliance, self._complex_love = \
            effective_rigid, complex_comp, complex_love

        if update_world_tides and not call_from_world:
            self.world.update_global_tides()

        return effective_rigid, complex_comp, complex_love

    def update_strength(self) -> Tuple[np.ndarray, np.ndarray]:
        """ Calculate viscosity and shear modulus and update the layer's state with new values.

        Calls self.viscosity_func and self.static_shear_modulus for pre-partial melting values.
        If self.partial_melt is present it will call upon it to calculate partial melted values.

        If pressure is set in the layer it will use it in the calculations (assuming that the viscosity and shear
        functions actually depend upon pressure). Otherwise, a pressure of zero will be assumed.

        Returns
        -------
        viscosity : np.ndarray
            Viscosity of the layer in [Pa s]
        shear_modulus : np.ndarray
            Shear Modulus of the layer in [Pa]
        """

        if self.pressure is None:
            pressure = np.zeros_like(self.temperature)
        else:
            pressure = self.pressure

        # Before Partial Melting is considered
        premelt_shear = self.static_shear_modulus
        premelt_visco = self.viscosity_func(self.temperature, pressure, *self.viscosity_inputs)
        liquid_viscosity = self.viscosity_liq_func(self.temperature, pressure, *self.viscosity_liq_inputs)

        # Partial Melt
        if self.partial_melt is not None:
            viscosity, shear_modulus = self.partial_melt.calculate(premelt_visco, premelt_shear, liquid_viscosity)
        else:
            viscosity, shear_modulus = premelt_visco, premelt_shear

        # Set state variables
        self._viscosity = viscosity
        self._shear_modulus = shear_modulus

        return viscosity, shear_modulus

    def calc_temperature_derivative(self) -> np.ndarray:
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

        for heat_source in self.heat_sources:
            if heat_source is None:
                raise ParameterMissingError
            total_heating += heat_source()

        total_cooling = self.heat_flux * self.surface_area_outer

        self._diff_temperature = (total_heating - total_cooling) / \
                                 (self.mass * self.specific_heat * (self.stefan + 1.) * self.temp_ratio)

        return self.diff_temperature

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

    # Class Properties
    @property
    def time(self):
        return self.world.time

    @time.setter
    def time(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_modes(self) -> Tuple[np.ndarray]:
        return self.world.tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value: Tuple[np.ndarray]):
        raise ImproperAttributeHandling('Tidal modes should be set at the world, not layer, level.')

    @property
    def tidal_freqs(self) -> Tuple[np.ndarray]:
        return self.world.tidal_freqs

    @tidal_freqs.setter
    def tidal_freqs(self, value: Tuple[np.ndarray]):
        raise ImproperAttributeHandling('Tidal frequencies should be set at the world, not layer, level.')

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
        if self.tidal_modes is not None:
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

    @property
    def melt_fraction(self) -> np.ndarray:

        if self._melt_fraction is None:
            raise AttributeNotSetError
        else:
            return self._melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        raise ImproperAttributeHandling('Melt fraction should be set by changing self.temperature or using by using'
                                        'the self.set_meltfraction method.')

    @property
    def viscosity(self) -> np.ndarray:

        if self._viscosity is None:
            raise AttributeNotSetError
        else:
            return self._viscosity

    @viscosity.setter
    def viscosity(self, value):

        if debug_mode:
            if type(value) not in floatarray_like:
                raise IncorrectAttributeType
            if type(value) == np.ndarray:
                if np.any(value > 1.0e30) or np.any(value < 1.0e-10):
                    raise UnusualRealValueError
            else:
                if value > 1.0e30 or value < 1.0e-10:
                    raise UnusualRealValueError

        self._viscosity = value_cleanup(value)

        # Update other models
        force_cooling = False
        force_tides = False
        if self.temperature is not None:
            force_cooling = True
            if self.tidal_modes is not None:
                force_tides = True
        self.update_tides(force_calculation=force_tides)
        self.update_cooling(force_calculation=force_cooling)

    @property
    def shear_modulus(self) -> np.ndarray:

        if self._shear_modulus is None:
            raise AttributeNotSetError
        else:
            return self._shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):

        if debug_mode:
            if type(value) not in floatarray_like:
                raise IncorrectAttributeType
            if type(value) == np.ndarray:
                if np.any(value > 1.0e20) or np.any(value < 1.0e-10):
                    raise UnusualRealValueError
            else:
                if value > 1.0e20 or value < 1.0e-10:
                    raise UnusualRealValueError

        self._shear_modulus = value_cleanup(value)

        # Update other models
        force_cooling = False
        force_tides = False
        if self.temperature is not None:
            force_cooling = True
            if self.tidal_modes is not None:
                force_tides = True
        self.update_tides(force_calculation=force_tides)
        self.update_cooling(force_calculation=force_cooling)

    # Attributes that will be updated whenever temperature is updated:
    @property
    def heat_flux(self):
        return self._heat_flux

    @heat_flux.setter
    def heat_flux(self, value):
        raise ImproperAttributeHandling('Heat flux is calculated whenever self.temperature is changed')

    @property
    def rayleigh(self):
        return self._rayleigh

    @rayleigh.setter
    def rayleigh(self, value):
        raise ImproperAttributeHandling('Rayleigh number is calculated whenever self.temperature is changed')

    @property
    def nusselt(self):
        return self._nusselt

    @nusselt.setter
    def nusselt(self, value):
        raise ImproperAttributeHandling('Nusselt number is calculated whenever self.temperature is changed')

    @property
    def blt(self):
        return self._blt

    @blt.setter
    def blt(self, value):
        raise ImproperAttributeHandling(
                'The thermal boundary layer thickness is calculated whenever self.temperature is changed')

    @property
    def complex_compliance(self) -> Tuple[np.ndarray]:
        return self._complex_compliance

    @complex_compliance.setter
    def complex_compliance(self, value):
        raise ImproperAttributeHandling

    @property
    def complex_love(self):
        return self._complex_love

    @complex_love.setter
    def complex_love(self, value) -> Tuple[np.ndarray]:
        raise ImproperAttributeHandling

    @property
    def effective_rigidity(self):
        return self._effective_rigidity

    @effective_rigidity.setter
    def effective_rigidity(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating(self):
        return self._tidal_heating

    @tidal_heating.setter
    def tidal_heating(self, value):

        if type(value) != np.ndarray:
            value = np.asarray(value)

        self._tidal_heating = value

    @property
    def radiogenic_heating(self):
        return self._radiogenic_heating

    @radiogenic_heating.setter
    def radiogenic_heating(self, value):

        if type(value) != np.ndarray:
            value = np.asarray(value)

        self._radiogenic_heating = value

    # Derivative Properties
    @property
    def diff_temperature(self):
        return self._diff_temperature

    @diff_temperature.setter
    def diff_temperature(self, value):
        raise ImproperAttributeHandling

    # Aliased Attributes
    @property
    def shear(self) -> np.ndarray:
        return self.shear_modulus

    @shear.setter
    def shear(self, value):
        self.shear_modulus = value

    @property
    def boundary_layer_thickness(self):
        return self.blt

    @boundary_layer_thickness.setter
    def boundary_layer_thickness(self, value):
        self.blt = value

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
