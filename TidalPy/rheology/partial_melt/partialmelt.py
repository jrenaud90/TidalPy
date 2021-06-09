from typing import Tuple, TYPE_CHECKING

import numpy as np

from . import known_model_live_args, known_model_const_args, known_models
from .defaults import partial_melt_defaults
from ... import log
from ...exceptions import BadValueError, InitiatedPropertyChangeError, IncorrectMethodToSetStateProperty, \
    OuterscopePropertySetError
from ...utilities.classes.model import LayerModelHolder
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray

if TYPE_CHECKING:
    from ..rheology import Rheology
    from ...structures.layers import PhysicalLayerType

@njit()
def calculate_melt_fraction(temperature: float, solidus: float, liquidus: float) -> float:
    """ Calculates the partial melt volume fraction based on the material's solidus and liquidus - NonArray Only

    Parameters
    ----------
    temperature : float
        Temperature of the layer or material [K]
    solidus : float
        Solidus temperature of the material [K]
    liquidus : float
        Liquidus temperature of the material [K]

    Returns
    -------
    partial_melt_volume_frac : float
        Volumetric Melt Fraction [m3 m-3]
    """

    if solidus >= liquidus:
        raise BadValueError('Solidus temperature can not be larger to equal to the liquidus temperature.')

    partial_melt_volume_frac = (temperature - solidus) / (liquidus - solidus)

    # Check for over/under-shoots
    if partial_melt_volume_frac < 0.:
        partial_melt_volume_frac = 0.
    elif partial_melt_volume_frac > 1.:
        partial_melt_volume_frac = 1.

    return partial_melt_volume_frac

@njit()
def calculate_melt_fraction_array(temperature: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the partial melt volume fraction based on the material's solidus and liquidus - Arrays Only

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

@njit()
def calculate_temperature_frommelt(melt_frac: float, solidus: float, liquidus: float) -> float:
    """ Calculates the temperature from the volumetric melt fraction - NonArray Only

    Parameters
    ----------
    melt_frac : float
        Volumetric Melt Fraction [m3 m-3]
    solidus : float
        Material/Layer solidus temperature [K]
    liquidus : float
        Material/Layer liquidus temperature [K]

    Returns
    -------
    temp_at_melt : float
        Temperature at melt fraction [K]

    """

    # Check for over/under-shoots
    if melt_frac < 0.:
        melt_frac = 0.
    elif melt_frac > 1.:
        melt_frac = 1.

    temp_at_melt = melt_frac * (liquidus - solidus) + solidus

    return temp_at_melt

@njit()
def calculate_temperature_frommelt_array(melt_frac: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the temperature from the volumetric melt fraction - Arrays Only

    Parameters
    ----------
    melt_frac : np.ndarray
        Volumetric Melt Fraction [m3 m-3]
    solidus : float
        Material/Layer solidus temperature [K]
    liquidus : float
        Material/Layer liquidus temperature [K]

    Returns
    -------
    temp_at_melt : np.ndarray
        Temperature at melt fraction [K]

    """

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    temp_at_melt = melt_frac * (liquidus - solidus) + solidus

    return temp_at_melt


class PartialMelt(LayerModelHolder):

    """ PartialMelt

    Partial melting provides a further temperature dependence to both viscosity and shear modulus. Depending upon the
        specific partial-melt model, the higher the temperature, the higher the partial melt fraction, and the lower
        the viscosity and shear modulus of the material.

    See Also
    --------
    TidalPy.utilities.methods.model.LayerModelHolder
    TidalPy.rheology.Rheology
    """

    default_config = partial_melt_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('rheology', 'partial_melting')

    def __init__(self, layer: 'PhysicalLayerType', rheology_class: 'Rheology', model_name: str = None,
                 store_config_in_layer: bool = True, initialize: bool = True):
        """ Constructor for ComplexCompliance class

        Parameters
        ----------
        layer : PhysicalLayerType
            The layer instance which the complex compliance should perform calculations on.
        rheology_class : Rheology
            Rheology class instance where the complex compliance model is stored.
        model_name : str = None
            The user-provided complex compliance name.
        store_config_in_layer: bool = True
            Flag that determines if the final complex compliance model's configuration dictionary should be copied
                into the `layer.config` dictionary.
        initialize : bool = True
            Determines if initial reinit should be performed on the model (loading in data from its `self.config`).
        """

        super().__init__(layer, model_name, store_config_in_layer, initialize=False)

        # Initialized properties
        self._rheology_class = rheology_class

        # Configuration properties
        self.solidus = None
        self.liquidus = None
        self.use_partial_melt = None

        # State properties
        self._melt_fraction = None
        self._postmelt_viscosity = None
        self._postmelt_shear_modulus = None
        self._postmelt_compliance = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Reinit method for the PartialMelting class

        Model will look at the user-provided configurations and pull out model information including constants

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this method has been called (additional steps may be
                preformed during the first reinit call).
        """

        # Load in configurations
        self.solidus = self.config.get('solidus', None)
        self.liquidus = self.config.get('liquidus', None)
        if self.model == 'off':
            self.use_partial_melt = False
        else:
            self.use_partial_melt = True

        if self.solidus is None or self.liquidus is None:
            log.debug(f'Solidus or Liquidus not provided in configurations for {self}.')

        super().reinit(initial_init)

    def clear_state(self):
        """ Clears the current state of the class without destroying data set by initialization.
        """

        super().clear_state()

        self._melt_fraction = None
        self._postmelt_viscosity = None
        self._postmelt_shear_modulus = None
        self._postmelt_compliance = None

    def _calculate(self) -> Tuple[FloatArray, FloatArray, FloatArray]:
        """ Wrapper for the partial melting function.

        First the partial melt will be updated based on the provided temperature. Then changes to the shear modulus
            and viscosity are updated.

        Returns
        -------
        melt_fraction : FloatArray
            Volumetric Melt Fraction [m3 m-3]
        postmelt_viscosity : FloatArray
            Post melting effective viscosity [Pa s]
        postmelt_shear_modulus : FloatArray
            Post melting shear modulus [Pa]
        """

        if self.use_partial_melt:
            if type(self.temperature) == np.ndarray:
                melt_fraction = calculate_melt_fraction_array(self.temperature, self.solidus, self.liquidus)
                postmelt_viscosity, postmelt_shear_modulus = \
                    self.func_array(melt_fraction, *self.live_inputs, *self.inputs)
            else:
                melt_fraction = calculate_melt_fraction(self.temperature, self.solidus, self.liquidus)
                postmelt_viscosity, postmelt_shear_modulus = \
                    self.func(melt_fraction, *self.live_inputs, *self.inputs)
        else:
            melt_fraction = np.zeros_like(self.temperature)
            postmelt_viscosity, postmelt_shear_modulus = self.premelt_viscosity, self.premelt_shear

        self._melt_fraction = melt_fraction
        self._postmelt_viscosity = postmelt_viscosity
        self._postmelt_shear_modulus = postmelt_shear_modulus
        self._postmelt_compliance = 1. / postmelt_shear_modulus

        return melt_fraction, postmelt_viscosity, postmelt_shear_modulus

    # Wrappers for user convenience
    def calculate_melt_fraction(self, temperature: FloatArray) -> FloatArray:

        if type(temperature) is np.ndarray:
            return calculate_melt_fraction_array(temperature, self.solidus, self.liquidus)
        else:
            return calculate_melt_fraction(temperature, self.solidus, self.liquidus)

    def calculate_temperature_frommelt(self, melt_fraction: FloatArray) -> FloatArray:

        if type(melt_fraction) is np.ndarray:
            return calculate_temperature_frommelt_array(melt_fraction, self.solidus, self.liquidus)
        else:
            return calculate_temperature_frommelt(melt_fraction, self.solidus, self.liquidus)


    # # Initialized properties
    @property
    def rheology_class(self) -> 'Rheology':
        """ The rheology class instance where the complex compliance model is stored """
        return self._rheology_class

    @rheology_class.setter
    def rheology_class(self, value):
        raise InitiatedPropertyChangeError


    # # State properties
    @property
    def postmelt_viscosity(self) -> FloatArray:
        """ Material's solid viscosity after the effects of partial melting have been applied """
        return self._postmelt_viscosity

    @postmelt_viscosity.setter
    def postmelt_viscosity(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def postmelt_shear_modulus(self) -> FloatArray:
        """ Material's shear modulus (rigidity) after the effects of partial melting have been applied """
        return self._postmelt_shear_modulus

    @postmelt_shear_modulus.setter
    def postmelt_shear_modulus(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def postmelt_compliance(self) -> FloatArray:
        """ Material's compliance (inverse of rigidity) after the effects of partial melting have been applied """
        return self._postmelt_compliance

    @postmelt_compliance.setter
    def postmelt_compliance(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def melt_fraction(self) -> FloatArray:
        """ Layer's volumetric melt fraction """
        return self._melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        raise IncorrectMethodToSetStateProperty


    # # Outer-scope properties
    @property
    def temperature(self):
        """ Outer-scope wrapper for `layer.temperature` """
        return self.layer.temperature

    @temperature.setter
    def temperature(self, value):
        raise OuterscopePropertySetError

    @property
    def premelt_viscosity(self):
        """ Outer-scope wrapper for `layer.premelt_viscosity` """
        return self.rheology_class.premelt_viscosity

    @premelt_viscosity.setter
    def premelt_viscosity(self, value):
        raise OuterscopePropertySetError

    @property
    def premelt_shear(self):
        """ Outer-scope wrapper for `layer.premelt_shear` """
        return self.rheology_class.premelt_shear

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise OuterscopePropertySetError

    @property
    def liquid_viscosity(self):
        """ Outer-scope wrapper for `layer.liquid_viscosity` """
        return self.rheology_class.liquid_viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        raise OuterscopePropertySetError
