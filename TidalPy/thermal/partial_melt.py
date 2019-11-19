import copy

import numpy as np

from TidalPy.rheology.partialmelt import melting_models
from .defaults import partial_melter_param_defaults
from .. import debug_mode
from ..types import FloatArray
from ..exceptions import BadAttributeValueError
from ..performance import njit
from ..utilities.model import LayerModel, ModelSearcher


# TODO: Move liquid shear into a material config?

# Partial Melt Model Finder
find_partial_melter = ModelSearcher(melting_models, partial_melter_param_defaults)


# Helper Functions
@njit
def calc_melt_fraction(temperature: FloatArray, solidus: float, liquidus: float) -> FloatArray:
    """ Calculates the volumetric melt fraction

    Parameters
    ----------
    temperature : FloatArray
        Material/Layer temperature [K]
    solidus : float
        Material/Layer solidus temperature [K]
    liquidus : float
        Material/Layer liquidus temperature [K]

    Returns
    -------
    melt_frac : FloatArray
        Volumetric Melt Fraction [m3 m-3]
    """

    melt_frac = (temperature - solidus) / (liquidus - solidus)

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    return melt_frac


@njit
def temp_from_melt(melt_frac: FloatArray, solidus: float, liquidus: float) -> FloatArray:
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


class PartialMelt(LayerModel):
    default_config = copy.deepcopy(partial_melter_param_defaults)
    config_key = 'partial_melt'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=find_partial_melter, call_reinit=True)

        self.temp_from_melt = temp_from_melt
        self.use_partial_melt = True
        self.solidus = None
        self.liquidus = None

        if self.config['model'] == 'off':
            self.use_partial_melt = False
            self.calculate = self._calculate_off

        if self.use_partial_melt:
            self.solidus = self.config['solidus']
            self.liquidus = self.config['liquidus']

    def calc_melt_fraction(self):
        """ Wrapper for calc_melt_fraction """

        if self.use_partial_melt:
            melt_fraction = calc_melt_fraction(self.layer.temperature, self.solidus, self.liquidus)
        else:
            melt_fraction = np.zeros_like(self.layer.temperature)

        if debug_mode:
            if np.any(melt_fraction > 1.) or np.any(melt_fraction < 0.):
                raise BadAttributeValueError

        return melt_fraction

    def _calculate(self, premelt_viscosity: np.ndarray, premelt_shear: np.ndarray, liquid_viscosity: np.ndarray):
        """ Wrapper for the partial melting function """

        return self.func(self.layer.temperature, self.layer.melt_fraction, premelt_viscosity, premelt_shear,
                         liquid_viscosity, *self.inputs)

    @staticmethod
    def _calculate_off(premelt_viscosity: np.ndarray, premelt_shear: np.ndarray, liquid_viscosity: np.ndarray):

        return premelt_viscosity, premelt_shear
