from __future__ import annotations

from TidalPy.radiogenics.radiogenics import Radiogenics
from TidalPy.structures.physical import PhysicalObjSpherical, ImproperAttributeChanged
from TidalPy.thermal.cooling import Cooling
from TidalPy.thermal.partial_melt import PartialMelt
from TidalPy.types import floatarray_like
from TidalPy.exceptions import (AttributeNotSetError, IncorrectAttributeType, UnusualRealValueError,
                                ImproperAttributeHandling,
                                BadAttributeValueError, UnknownTidalPyConfigValue, ReinitError)
from ..thermal import find_viscosity, calc_melt_fraction
from ..rheology.rheology import Rheology
import numpy as np
from typing import Dict, List
from scipy.constants import G
from .. import log, debug_mode
from typing import Tuple
import burnman
from ..utilities.numpy_help import find_nearest
from ..bm.conversion import burnman_property_name_conversion, burnman_property_value_conversion
from ..configurations import burnman_interpolation_method, burnman_interpolation_N
from TidalPy.utilities.dict_tools import nested_get
from typing import Union
from .defaults import layer_defaults
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .worlds import BurnManWorld

class ThermalLayer(PhysicalObjSpherical):

    """ Layer Class: Tidal Planets are made up of various Layers of different compositions.

        Any functionality that requires information about material properties should be called from the layer class.
    """

    default_config = layer_defaults


    def __init__(self, layer_name: str, world: BurnManWorld, burnman_layer: burnman.Layer, layer_config: dict):

        # Load layer defaults
        self.type = layer_config['type']
        self.default_config = self.default_config[self.type]

        # Setup Physical Layer Geometry
        super().__init__(layer_config, automate=True)

        self.name = layer_name
        self.world = world

        # Information about the other layers (mostly for gravity and cooling calculations)
        self.layer_below = None             # type: ThermalLayer or None
        self.layer_above = None             # type: ThermalLayer or None

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
        self.gravity_surf = self.bm_layer.gravity[-1]

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
        self.thermal_diffusivity = None   # type: float or None
        self.thermal_conductivity = None  # type: float or None
        self.heat_fusion = None           # type: float or None
        self.temp_ratio = None            # type: float or None
        self.stefan = None                # type: float or None

        # State Variables
        self._temperature = None          # type: np.ndarray or None
        self._melt_fraction = None        # type: np.ndarray or None
        self._viscosity = None            # type: np.ndarray or None
        self._shear_modulus = None        # type: np.ndarray or None
        self._heat_flux = None            # type: np.ndarray or None
        self._rayleigh = None             # type: np.ndarray or None
        self._nusselt = None              # type: np.ndarray or None
        self._blt = None                  # type: np.ndarray or None
        self._temperature_lower = None    # type: np.ndarray or None

        # Model Holders
        self.cooling = None               # type: Cooling or None
        self.partial_melt = None          # type: PartialMelt or None
        self.radiogenics = None           # type: Radiogenics or None
        self.heat_sources = list()

        # Function Holders
        self.viscosity_func, self.viscosity_inputs, self.viscosity_live_inputs = None, None, None
        self.viscosity_liq_func, self.viscosity_liq_inputs,self.viscosity_liq_live_inputs = None, None, None

        # Call to init
        self.init()

    def _build_material_property_interpolation(self):
        """ Interpolates material properties based on a fixed pressure and a suggested temperature range.

            This can take some time to precompute
        """

        if self.pressure is None:
            raise AttributeNotSetError

        interp_properties = ['bulk_modulus', 'thermal_expansion', 'specific_heat']
        bm_properties = [burnman_property_name_conversion[interp_prop] for interp_prop in interp_properties]

        # We will use the fixed pressure to calculate the various parameters at all the temperature ranges
        #     These results will then be used in an interpolation for whenever self.temperature changes
        pressures = self.pressure * np.ones_like(self.interp_temperature_range)
        property_results = None
        try:
            property_results = self.bm_material.evaluate(bm_properties, pressures, self.interp_temperature_range)
        except TypeError as e:
            # FIXME: There is some bug in bug in Burnman when it trys to send a warning. If that is what got caught, ignore it for now.
            if '<lambda>() got an unexpected keyword argument' in e.args[0]:
                pass
            else:
                raise e

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

    def init(self):

        # Material properties that might have been affected by new configuration files
        self.static_shear_modulus = self.config['shear_modulus']
        self.heat_fusion = self.config['heat_fusion']
        self.thermal_conductivity = self.config['thermal_conductivity']
        self.temp_ratio = self.config['boundary_temperature_ratio']
        self.stefan = self.config['stefan']

        # Setup viscosity functions
        solid_visco_model = nested_get(self.config, ['rheology', 'solid_viscosity', 'model'], default=None)
        liquid_visco_model = nested_get(self.config, ['rheology', 'liquid_viscosity', 'model'], default=None)

        self.viscosity_func, self.viscosity_inputs, self.viscosity_live_inputs = \
            find_viscosity(solid_visco_model, default_key=[self.type, 'solid_viscosity'])
        self.viscosity_liq_func, self.viscosity_liq_inputs, self.viscosity_liq_live_inputs = \
            find_viscosity(liquid_visco_model, default_key=[self.type, 'liquid_viscosity'])

        # Setup Partial Melting
        self.partial_melt = PartialMelt(self)

        # Setup Cooling
        self.cooling = Cooling(self)

        # Setup Radiogenics
        self.radiogenics = Radiogenics(self)
        self.heat_sources.append(self.radiogenics)

        super().init()

    def reinit(self):

        # Check for changes to config that may break the planet
        for layer_name, layer_dict in self.config['layers'].items():
            old_layer_dict = self._old_config[layer_name]
            for critical_attribute in ['material', 'material_source', 'type', 'thickness', 'radius', 'slices']:
                if layer_dict[critical_attribute] != old_layer_dict[critical_attribute]:
                    raise ReinitError

        super().reinit()
        self.set_geometry(self.radius, self.mass, self.thickness)

    def set_geometry(self, radius: float, mass: float, thickness: float = None):

        super().set_geometry(radius, mass, thickness)

    def set_strength(self, viscosity: np.ndarray = None, shear_modulus: np.ndarray = None):
        """ Sets the strength via viscosity and shear modulus of the layer.

        Use this instead of the viscosity.setter or shear_modulus.setter when you are changing both at the same time.
        It will use less calls.
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
            if type(viscosity) != np.ndarray:
                viscosity = np.asarray(viscosity)
            self._viscosity = viscosity
        if shear_modulus is not None:
            if type(shear_modulus) != np.ndarray:
                shear_modulus = np.asarray(shear_modulus)
            self._shear_modulus = shear_modulus

        # Calculate cooling
        self._heat_flux, self._blt, self._rayleigh, self._nusselt = self.cooling.calculate()


    def update_strength(self) -> Tuple[np.ndarray, np.ndarray]:
        """ Sets and returns the strength (viscosity and shear) of the layer

        Requires self.viscosity_func and self.shear_modulus
        Optionally it will call self.partial_melt_func if provided

        :return: <Tuple[FloatOrArray]> (Viscosity, Shear Modulus)
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
        viscosity, shear_modulus = self.partial_melt.calculate(premelt_visco, premelt_shear, liquid_viscosity)

        # Set state variables
        self._viscosity = viscosity
        self._shear_modulus = shear_modulus

        return viscosity, shear_modulus

    def find_diff_temperature(self):

        # Initialize total heating based on if there is a layer warming this one from below
        if self.layer_below is None:
            total_heating = np.zeros_like(self.temperature)
        else:
            total_heating = self.layer_below.heat_flux * self.surface_area_inner

        for heat_source in self.heat_sources:
            if 'calculate' in heat_source.__dict__:
                total_heating += heat_source.calculate()
            else:
                total_heating += heat_source()

        cooling_flux, cooling_thickness, rayleigh, nusselt = self.cooling.calculate()
        total_cooling = cooling_flux * self.surface_area_outer

        return (total_heating - total_cooling) / self.energy_per_therm

    # Class Properties
    @property
    def temperature(self) -> np.ndarray:

        if self._temperature is None:
            raise AttributeNotSetError
        else:
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
        if type(value) != np.ndarray:
            value = np.asarray(value)
        self._temperature = value

        # Update shear and viscosity
        self.update_strength()

        # Set new melt fraction
        self._melt_fraction = self.partial_melt.calc_melt_fraction()
        # Set new material properties based on BurnMan Interpolation
        for prop_name, prop_data in self._interp_prop_data_lookup.items():
            setattr(self, prop_name, np.interp(self._temperature, self.interp_temperature_range, prop_data))

        # Set new dependent material properties
        self.thermal_diffusivity = self.thermal_conductivity / (self.density * self.specific_heat)
        self.energy_per_therm = self.mass * self.specific_heat * (self.stefan + 1.) * self.temp_ratio

        # Calculate cooling
        self._heat_flux, self._blt, self._rayleigh, self._nusselt = self.cooling.calculate()

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
        raise ImproperAttributeHandling('The thermal boundary layer thickness is calculated whenever self.temperature is changed')

    # TODO: For now lower temperature just returns self.temperature. This is used in cooling calculations
    @property
    def temperature_lower(self):
        return self.temperature

    @temperature_lower.setter
    def temperature_lower(self, value):
        raise ImproperAttributeHandling('temperature_lower is calculated by setting self.temperature.')

    @property
    def temperature_surf(self):
        if self.layer_above is not None:
            return self.layer_above.temperature_lower
        else:
            return self.config['surface_temperature']

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
        if type(value) != np.ndarray:
            value = np.asarray(value)
        self._viscosity = value

        # Calculate cooling
        self._heat_flux, self._blt, self._rayleigh, self._nusselt = self.cooling.calculate()

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
        if type(value) != np.ndarray:
            value = np.asarray(value)
        self._shear_modulus = value

        # Calculate cooling
        self._heat_flux, self._blt, self._rayleigh, self._nusselt = self.cooling.calculate()

    # Aliased Attributes
    @property
    def shear(self) -> np.ndarray:
        return self.shear_modulus

    @shear.setter
    def shear(self, value):
        self.shear_modulus = value


class TidalLayer(ThermalLayer):

    def __init__(self, layer_name: str, world, burnman_layer: burnman.Layer, layer_config: dict):

        super().__init__(layer_name, world, burnman_layer, layer_config)


        # Constants
        self.tidal_scale = 1.
        if self.config['use_tvf']:
            self.tidal_scale = self.volume / self.world.volume

        # Models
        self.rheology = None             # type: Rheology or None

        # State variables
        self._complex_compliance = None  # type: np.ndarray or None
        self._complex_love = None        # type: np.ndarray or None
        self._effective_rigidity = None  # type: np.ndarray or None
        self._tidal_heating = None       # type: np.ndarray or None

    def init(self):

        super().init()
        self.rheology = Rheology(self)

    def set_strength(self, viscosity: np.ndarray = None, shear_modulus: np.ndarray = None):

        super().set_strength(viscosity, shear_modulus)
        self.world.update_tides()

    def find_diff_temperature(self):

        diff_temp_notides = super().find_diff_temperature()

        # Now include tides
        diff_temp_dt = {rheo_name: diff_temp_notides + tidal_heating / self.energy_per_therm for
                        rheo_name, tidal_heating in self.tidal_heating.items()}

        return diff_temp_dt

    @property
    def tidal_modes(self) -> Tuple[np.ndarray]:
        return self.world.tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value: Tuple[np.ndarray]):
        raise ImproperAttributeHandling('Tidal modes should be set at the world, not layer, level.')

    @property
    def complex_compliance(self) -> Dict[str, Tuple[np.ndarray]]:
        return self._complex_compliance

    @complex_compliance.setter
    def complex_compliance(self, value):
        raise ImproperAttributeHandling

    @property
    def complex_love(self):
        return self._complex_love

    @complex_love.setter
    def complex_love(self, value) -> Dict[str, Tuple[np.ndarray]]:
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
        raise ImproperAttributeHandling

    # Update previous class state variables to include new setter functionality
    @property
    def temperature(self) -> np.ndarray:
        return super().temperature

    @temperature.setter
    def temperature(self, value):
        super(TidalLayer, self.__class__).temperature.fset(self, value)

        self._effective_rigidity, self._complex_compliance, self._complex_love = \
            self.rheology.calculate()

        self.world.update_tides()

    @property
    def viscosity(self) -> np.ndarray:
        return super().viscosity

    @viscosity.setter
    def viscosity(self, value):

        super(TidalLayer, self.__class__).viscosity.fset(self, value)

        self._effective_rigidity, self._complex_compliance, self._complex_love = \
            self.rheology.calculate()

        self.world.update_tides()

    @property
    def shear_modulus(self) -> np.ndarray:
        return super().shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        super(TidalLayer, self.__class__).shear_modulus.fset(self, value)

        self._effective_rigidity, self._complex_compliance, self._complex_love = \
            self.rheology.calculate()

        self.world.update_tides()

# Helpers
def construct_layer(layer_name: str, world: BurnManWorld, burnman_layer: burnman.Layer, layer_config: dict):

    # Try to determine if the layer is tidal or not
    is_tidal = layer_config['is_tidal']

    if 'rheology' in layer_config:
        is_tidal = True

    if is_tidal:
        return TidalLayer(layer_name, world, burnman_layer, layer_config)
    else:
        return ThermalLayer(layer_name, world, burnman_layer, layer_config)