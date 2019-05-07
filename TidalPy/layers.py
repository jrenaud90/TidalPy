from .physical import PhysicalObjSpherical
from .types import floatarray_like
from .exceptions import AttributeNotSet, IncorrectAttributeType, UnusualRealValueError, ImproperAttributeHandling,\
    BadAttributeValueError
from . import debug_mode
from .thermal import melt_fraction, Strength
import numpy as np


class Layer(PhysicalObjSpherical):

    """ Layer Class: Tidal Planets are made up of various Layers of different compositions.

        Any functionality that requires information about material properties should be called from the layer class.
    """

    def __init__(self, bm_layer: burnman.Layer, layer_below: Layer, material_properties: dict):




        self.layer_below = layer_below
        self.mass_below = layer_below.mass

        # Pull out information from the already initialized burnman Layer
        self.bm_layer = bm_layer

        # State Variables
        self._temperature = None
        self._melt_fraction = None
        self._strength = None

        # Material Properties

        # Setup Viscosity and Shear Functions
        # TODO
        self.viscosity_func =
        self.shear_func =

        # Setup Partial Melting functions
        self.use_partial_melt = None
        self.solidus = None
        self.liquidus = None
        # TODO
        self.partial_melt_func =

    def calc_strength(self) -> Strength:
        """ Sets and returns the material strength (viscosity and shear modulus)

        requires the temperature of the layer to be set

        """

        pre_melt_visco = self.viscosity_func(self.temperature)
        pre_melt_shear = self.shear_func(self.temperature)

        if not self.use_partial_melt:
            return Strength(pre_melt_visco, pre_melt_shear)
        else:
            return Strength(*self.partial_melt_func(melt_fraction, pre_melt_visco, pre_melt_shear))

    def set_meltfraction(self) -> np.ndarray:
        """ Sets and returns the melt fraction of the layer

        requires self.solidus and self.liquidus to be set
        :return: <FloatOrArray> Volumetric Melt Fraction
        """

        if self.solidus is None or self.liquidus is None:
            raise AttributeNotSet

        self._melt_fraction = melt_fraction(self.temperature, self.solidus, self.liquidus)

        if debug_mode:
            if np.any(self._melt_fraction > 1.) or np.any(self._melt_fraction < 1.):
                raise BadAttributeValueError

        return self.melt_fraction

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
        if self.use_partial_melt:
            self.set_meltfraction()

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


