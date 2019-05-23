from TidalPy.structures.physical import PhysicalObjSpherical, ImproperAttributeChanged
from TidalPy.thermal.cooling import Cooling
from TidalPy.thermal.partial_melt import PartialMelt
from TidalPy.types import floatarray_like
from TidalPy.exceptions import (AttributeNotSetError, IncorrectAttributeType, UnusualRealValueError,
                                ImproperAttributeHandling,
                                BadAttributeValueError, UnknownTidalPyConfigValue)
from TidalPy import debug_mode
from TidalPy.thermal import find_viscosity, calc_melt_fraction
import numpy as np
from scipy.constants import G
from typing import Tuple
import burnman
from ..utilities.numpy_help import find_nearest
from ..bm.conversion import burnman_property_name_conversion, burnman_property_value_conversion
from ..configurations import burnman_interpolation_method
from TidalPy.utilities.dict_tools import nested_get

class ThermalLayer(PhysicalObjSpherical):

    """ Layer Class: Tidal Planets are made up of various Layers of different compositions.

        Any functionality that requires information about material properties should be called from the layer class.
    """

    def __init__(self, layer_name: str, burnman_layer: burnman.Layer, mass_below, layer_config: dict):

        super().__init__(layer_config, automate=True)

        self.name = layer_name

        # Information about the other layers (mostly for gravity and cooling calculations)
        self.layer_below = None             # type: ThermalLayer
        self.layer_above = None             # type: ThermalLayer
        self.mass_below = mass_below

        # Pull out information from the already initialized burnman Layer
        self.bm_layer = burnman_layer
        self.bm_material = burnman_layer.material

        # Other Attributes - these should not change on a reinit
        self.pressure_upper = self.bm_layer.pressures[-1]
        self.pressure_lower = self.bm_layer.pressures[0]
        self.density_upper = self.bm_layer.density[-1]
        self.density_lower = self.bm_layer.density[0]
        self.bm_mid_index = find_nearest()

        # Attributes calculated by BurnMan but set at a specific spot
        if burnman_interpolation_method == 'mid':
            self.pressure = self.bm_layer.pressures[self.bm_mid_index]
            self.density = self.bm_layer.density[self.bm_mid_index]
            self.interp_func = lambda array: array[self.bm_mid_index]
        elif burnman_interpolation_method == 'avg':
            self.pressure = np.average(self.bm_layer.pressures)
            self.density = np.average(self.bm_layer.density)
            self.interp_func = np.average
        elif burnman_interpolation_method == 'median':
            self.pressure = np.median(self.bm_layer.pressures)
            self.density = np.median(self.bm_layer.density)
            self.interp_func = np.median
        else:
            raise UnknownTidalPyConfigValue

        # Setup geometry based on the BurnMan results
        self.gravity = None
        bm_radius = np.max(self.bm_layer.radii)
        bm_thickness = bm_radius - np.min(self.bm_layer.radii)
        bm_mass = self.bm_layer.mass
        self.set_geometry(radius=bm_radius, mass=bm_mass, thickness=bm_thickness)
        self.bm_mid_index = find_nearest(self.bm_layer.radii, self.radius - self.thickness/2)

        # Setup Physical Layer Geometry
        self.type = self.config['type']
        self.slices = self.config['slices']
        self.material = self.config['material']

        # Material Properties
        self.bulk_modulus =
        self.thermal_expansion =
        self.specific_heat =

        # Material properties that might have been affected by new configuration file (BurnMan does not calculate these)
        self.static_shear_modulus = None
        self.thermal_diffusivity = None
        self.thermal_conductivity = None
        self.heat_fusion = None

        # State Variables
        self._temperature = None
        self._melt_fraction = None
        self._viscosity = None
        self._shear_modulus = None

        # Model Holders
        self.cooling = None       # type: Cooling
        self.partial_melt = None  # type: PartialMelt

        # Function Holders
        self.viscosity_func, self.viscosity_inputs = None, None
        self.viscosity_liq_func, self.viscosity_liq_inputs = None, None

        # Call to init
        self.init()

    def _build_material_property_interpolation(self, property_name: str, temperature_range):
        """ Interpolates material properties based on a fixed pressure and a suggested temperature range.

            This can take some time to precompute
        """

        if self.pressure is None:
            raise AttributeNotSetError

        interp_properties = ['bulk_modulus', 'thermal_expansion', 'specific_heat']
        bm_properties = [burnman_property_name_conversion[interp_prop] for interp_prop in interp_properties]



        for interp_prop in ['bulk_modulus', 'thermal_expansion', 'specific_heat']:
            bm_name = burnman_property_name_conversion[interp_prop]

            conversion = 1.
            if interp_prop in burnman_property_value_conversion:
                if burnman_property_value_conversion[interp_prop] == 'molar':
                    conversion = 1. / self.bm_material.molar_mass
                else:
                    raise KeyError




    def init(self):

        # Interpolation Method
        self.interp_meth = self.config['interp_method']

        # Pull out information from the already initialized burnman Layer
        self.pressure = None
        bm_radius =
        bm_thickness =
        bm_mass =
        # TODO: For now we do not assume any pre-melt function for shear_modulus. Use BM material to set it. If we want to implement it in the future: use the same method used for viscosity
        self.shear_modulus_static =

        # Material properties that might have been affected by new configuration files
        self.static_shear_modulus = self.config['shear_modulus']
        self.heat_fusion = self.config['heat_fusion']
        self.thermal_diffusivity = self.thermal_conductivity / (self.density * self.specific_heat)

        # Setup Physical Layer Geometry
        self.type = self.config['type']


        # Material Properties
        #TODO
        self.use_bm_mat_props = self.config.get('use_burnman_material_properties', True)
        if self.use_bm_mat_props:
            self.config['convection']['thermal_conductivity'] =

        # Setup viscosity functions
        solid_visco_model = nested_get(self.config, ['rheology', 'solid_viscosity', 'model'], default=None)
        liquid_visco_model = nested_get(self.config, ['rheology', 'liquid_viscosity', 'model'], default=None)

        self.viscosity_func, self.viscosity_inputs = \
            find_viscosity(solid_visco_model, default_key=[self.type, 'solid_viscosity'])
        self.viscosity_liq_func, self.viscosity_liq_inputs = \
            find_viscosity(liquid_visco_model, default_key=[self.type, 'liquid_viscosity'])

        # Setup partial melting
        self.partial_melt = PartialMelt(self.type, self.config['partial_melt'])

        # Setup Cooling
        cooling_config = self.config['cooling']
        self.cooling = Cooling(self.type, self.thickness, cooling_config, self.gravity, self.density_bulk)
        self.calc_cooling = self.cooling.calculate

    def reinit(self):

        # Check for changes to config that may break the planet
        for critical_attribute in ['type', 'slices', 'material', 'thickness', 'radius']:
            if getattr(self, critical_attribute) != self.config[critical_attribute]:
                raise ImproperAttributeChanged

        super().reinit()
        self.set_geometry(self.radius, self.mass, self.thickness)

    def set_geometry(self, radius: float, mass: float, thickness: float = None):

        super().set_geometry(radius, mass, thickness)
        self.gravity = G * (self.mass + self.mass_below) / self.radius**2


    def set_meltfraction(self) -> np.ndarray:
        """ Sets and returns the melt fraction of the layer

        requires self.solidus and self.liquidus to be set
        :return: <FloatArray> Volumetric Melt Fraction
        """

        melt_fraction = self.partial_melt.calc_melt_fraction(self.temperature)
        if debug_mode:
            if np.any(melt_fraction > 1.) or np.any(melt_fraction < 0.):
                raise BadAttributeValueError

        self._melt_fraction = melt_fraction
        return melt_fraction

    def set_strength(self) -> Tuple[np.ndarray, np.ndarray]:
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
        premelt_shear = self.shear_modulus_static
        premelt_visco = self.viscosity_func(self.temperature, pressure, *self.viscosity_inputs)
        liquid_viscosity = self.viscosity_liq_func(self.temperature, pressure, *self.viscosity_liq_inputs)

        # Partial Melt
        viscosity, shear_modulus = self.partial_melt.calculate(self.temperature, self.melt_fraction,
                                                               premelt_visco, premelt_shear, liquid_viscosity)

        # Set state variables
        self._viscosity = viscosity
        self._shear_modulus = shear_modulus

        return viscosity, shear_modulus

    def find_cooling(self):
        """ Calculates cooling based on the states of this layer and the layer above it. """

        if debug_mode:
            if self.layer_above is None:
                raise AttributeNotSetError('In order to use a layer.find_cooling, the layer.layer_above must be set.')
        top_temp = self.layer_above.temperature

        return self.cooling.calculate(self.temperature, top_temp, self.viscosity, self.thickness)

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
        self.set_meltfraction()
        self.set_strength()

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
        raise ImproperAttributeHandling('Viscosity should be set by changing self.temperature or using by using'
                                        'the self.set_strength method.')

    @property
    def shear_modulus(self) -> np.ndarray:

        if self._shear_modulus is None:
            raise AttributeNotSetError
        else:
            return self._shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        raise ImproperAttributeHandling('Shear Modulus should be set by changing self.temperature or using by using'
                                        'the self.set_strength method.')

    # Aliased Attributes
    @property
    def shear(self) -> np.ndarray:
        return self.shear_modulus

    @shear.setter
    def shear(self, value):
        self.shear_modulus = value


class TidalLayer(ThermalLayer):
    pass