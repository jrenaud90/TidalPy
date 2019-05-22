from TidalPy.structures.physical import PhysicalObjSpherical
from TidalPy.thermal.cooling import Cooling
from TidalPy.thermal.partial_melt import PartialMelt
from TidalPy.types import floatarray_like
from TidalPy.exceptions import AttributeNotSet, IncorrectAttributeType, UnusualRealValueError, ImproperAttributeHandling,\
    BadAttributeValueError
from TidalPy import debug_mode
from TidalPy.thermal import find_viscosity, calc_melt_fraction
import numpy as np
from typing import Tuple

from TidalPy.utilities.dict_tools import nested_get


class ThermalLayer(PhysicalObjSpherical):

    """ Layer Class: Tidal Planets are made up of various Layers of different compositions.

        Any functionality that requires information about material properties should be called from the layer class.
    """

    def __init__(self, bm_layer: burnman.Layer, mass_below, layer_config: dict):

        super().__init__(layer_config, automate=True)

        if not isinstance(layer_below, (ThermalLayer, TidalLayer)):
            raise TypeError

        self.layer_below = None
        self.layer_above = None
        self.mass_below = mass_below

        # Pull out information from the already initialized burnman Layer
        self.bm_layer = bm_layer
        self.pressure = None
        bm_radius =
        bm_thickness =
        bm_mass =
        # TODO: For now we do not assume any pre-melt function for shear_modulus. Use BM material to set it. If we want to implement it in the future: use the same method used for viscosity
        self.shear_modulus_static =


        # Setup Physical Layer Geometry
        self.type = self.config['type']
        self.set_geometry(radius=bm_radius, thickness=bm_thickness, mass=bm_mass)

        # State Variables
        self._temperature = None
        self._melt_fraction = None
        self._viscosity = None
        self._shear_modulus = None

        # Material Properties
        #TODO
        self.use_bm_mat_props = self.config.get('use_burnman_material_properties', True)
        if self.use_bm_mat_props:
            self.config['convection']['thermal_conductivity'] =

        self.cooling = Cooling(self.type, self.thickness, cooling_config, self.gravity, self.density_bulk)

        # Setup viscosity functions
        solid_visco_model = nested_get(self.config, ['rheology', 'solid_viscosity', 'model'], default=None)
        liquid_visco_model = nested_get(self.config, ['rheology', 'liquid_viscosity', 'model'], default=None)

        self.viscosity_func, self.viscosity_inputs = \
            find_viscosity(solid_visco_model, default_key=[self.type, 'solid_viscosity'])
        self.viscosity_liq_func, self.viscosity_liq_inputs = \
            find_viscosity(liquid_visco_model, default_key=[self.type, 'liquid_viscosity'])

        # Setup partial melting
        self.partial_melt = PartialMelt(self.type, self.config['partial_melt'])

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

    def calc_rayleigh(self):
        """ Calculate Rayleigh Number """

    # Class Properties
    @property
    def temperature(self) -> np.ndarray:

        if self._temperature is None:
            raise AttributeNotSet
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
            raise AttributeNotSet
        else:
            return self._melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        raise ImproperAttributeHandling('Melt fraction should be set by changing self.temperature or using by using'
                                        'the self.set_meltfraction method.')

    @property
    def viscosity(self) -> np.ndarray:

        if self._viscosity is None:
            raise AttributeNotSet
        else:
            return self._viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise ImproperAttributeHandling('Viscosity should be set by changing self.temperature or using by using'
                                        'the self.set_strength method.')

    @property
    def shear_modulus(self) -> np.ndarray:

        if self._shear_modulus is None:
            raise AttributeNotSet
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