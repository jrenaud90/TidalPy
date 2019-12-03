from typing import TYPE_CHECKING

import numpy as np

from ...exceptions import BadValueError, ImproperAttributeHandling
from ...performance import njit
from ...utilities.model import LayerModelHolder
from . import known_model_live_args, known_model_const_args, known_models
from .defaults import partial_melt_defaults

if TYPE_CHECKING:
    from ...structures.layers import ThermalLayer
    from ..rheology import Rheology


@njit
def calculate_melt_fraction(temperature: np.ndarray, solidus: float, liquidus: float):
    """ Calculates the partial melt volume fraction based on the material's solidus and liquidus

    Parameters
    ----------
    temperature : np.ndarray
        Temperature of the layer or material [K]
    solidus : float
        Solidus temperature of the material [K]
    liquidus : float
        Liquidus temperature of the material [K]

    Returns
    -------
    partial_melt_volume_frac : np.ndarray
        Volumetric Melt Fraction [m3 m-3]
    """

    if solidus >= liquidus:
        raise BadValueError('Solidus temperature can not be larger to equal to the liquidus temperature.')

    partial_melt_volume_frac = (temperature - solidus) / (liquidus - solidus)

    # Check for over/under-shoots
    partial_melt_volume_frac[partial_melt_volume_frac < 0.] = 0.
    partial_melt_volume_frac[partial_melt_volume_frac > 1.] = 1.

    return partial_melt_volume_frac

@njit
def calculate_temperature_frommelt(melt_frac: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the temperature from the volumetric melt fraction

    Parameters
    ----------
    melt_frac : FloatArray
        Volumetric Melt Fraction [m3 m-3]
    solidus : float
        Material/Layer solidus temperature [K]
    liquidus : float
        Material/Layer liquidus temperature [K]

    Returns
    -------
    temp_at_melt : FloatArray
        Temperature at melt fraction [K]

    """

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    temp_at_melt = melt_frac * (liquidus - solidus) + solidus

    return temp_at_melt


class PartialMelt(LayerModelHolder):

    """ Partial Melting Class - Child of LayerModelHolder Class

    Partial melting provides a further temperature dependence to both viscosity and shear modulus. Depending upon the
        specific partial-melt model, the higher the temperature, the higher the partial melt fraction, and the lower
        the viscosity and shear modulus of the material.
    """

    default_config = partial_melt_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('rheology', 'partial_melting')

    def __init__(self, layer: ThermalLayer, rheology_class: Rheology, model_name: str = None,
                 store_config_in_layer: bool = True):

        super().__init__(layer, model_name, store_config_in_layer)

        self.rheology_class = rheology_class

        # Pull out specific information related to this module
        self.solidus = self.get_param('solidus')
        self.liquidus = self.get_param('liquidus')
        self.use_partial_melt = True
        if self.model == 'off':
            self.use_partial_melt = False

        # State attributes
        self._melt_fraction = None
        self._postmelt_viscosity = None
        self._postmelt_shear_modulus = None
        self._postmelt_compliance = None

    def _calculate(self):
        """ Wrapper for the partial melting function.

        First the partial melt will be updated based on the provided temperature. Then changes to the shear modulus
            and viscosity are updated.

        Returns
        -------
        melt_fraction : np.ndarray
            Volumetric Melt Fraction [m3 m-3]
        postmelt_viscosity : np.ndarray
            Post melting effective viscosity [Pa s]
        postmelt_shear_modulus : np.ndarray
            Post melting shear modulus [Pa]
        """

        if self.use_partial_melt:
            melt_fraction = calculate_melt_fraction(self.temperature, self.solidus, self.liquidus)
            postmelt_viscosity, postmelt_shear_modulus = self.func(melt_fraction, *self.live_inputs, *self.inputs)
        else:
            melt_fraction = np.zeros_like(self.temperature)
            postmelt_viscosity, postmelt_shear_modulus = self.premelt_viscosity, self.premelt_shear

        self._melt_fraction = melt_fraction
        self._postmelt_viscosity = postmelt_shear_modulus
        self._postmelt_shear_modulus = postmelt_shear_modulus
        self._postmelt_compliance = 1. / postmelt_shear_modulus

        return melt_fraction, postmelt_viscosity, postmelt_shear_modulus

    # Wrappers for user convenience
    def calculate_melt_fraction(self, temperature: np.ndarray):
        return calculate_melt_fraction(temperature, self.solidus, self.liquidus)

    def calculate_temperature_frommelt(self, melt_fraction: np.ndarray):
        return calculate_temperature_frommelt(melt_fraction, self.solidus, self.liquidus)


    # State properties
    @property
    def postmelt_viscosity(self):
        return self._postmelt_viscosity

    @postmelt_viscosity.setter
    def postmelt_viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def postmelt_shear_modulus(self):
        return self._postmelt_shear_modulus

    @postmelt_shear_modulus.setter
    def postmelt_shear_modulus(self, value):
        raise ImproperAttributeHandling

    @property
    def postmelt_compliance(self):
        return self._postmelt_compliance

    @postmelt_compliance.setter
    def postmelt_compliance(self, value):
        raise ImproperAttributeHandling

    @property
    def melt_fraction(self):
        return self._melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        raise ImproperAttributeHandling

    # Outerscope properties
    @property
    def temperature(self) -> np.ndarray:
        return self.layer.temperature

    @temperature.setter
    def temperature(self, value):
        raise ImproperAttributeHandling

    @property
    def premelt_viscosity(self) -> np.ndarray:
        return self.rheology_class.premelt_viscosity

    @premelt_viscosity.setter
    def premelt_viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def premelt_shear(self) -> np.ndarray:
        return self.rheology_class.premelt_shear

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise ImproperAttributeHandling

    @property
    def liquid_viscosity(self) -> np.ndarray:
        return self.rheology_class.liquid_viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        raise ImproperAttributeHandling
