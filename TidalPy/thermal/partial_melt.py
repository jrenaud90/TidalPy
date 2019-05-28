from numba import njit
import numpy as np

from TidalPy import debug_mode
from TidalPy.exceptions import BadAttributeValueError
from TidalPy.structures.layers import LayerType
from TidalPy.utilities.classes import ModelHolder, LayerModel
from ..utilities.search import ModelSearcher
from . import viscosity

# TODO: Move liquid shear into a material config?

# Partial Melt Model Finder
partial_melter_param_defaults = {
    'ice': {
        'model': 'off',
        'solidus': 270.,
        'liquidus': 273.15,
        # Liquid Shear is just something very small.
        'liquid_shear': 1.0e-5,
    },
    'rock': {
        'model': 'henning',
        'solidus': 1600.,
        'liquidus': 2000.,
        'liquid_shear': 1.0e-5,
        'fs_visc_power_slope': 27000.0,
        'fs_visc_power_phase': 1.0,
        'fs_shear_power_slope': 82000.0,
        'fs_shear_power_phase': 40.6,
        'crit_melt_frac': 0.5,
        'crit_melt_frac_width': 0.05,
        'hn_visc_slope_1': 13.5,
        'hn_visc_slope_2': 370.0,
        'hn_shear_param_1': 40000.0,
        'hn_shear_param_2': 25.0,
        'hn_shear_falloff_slope': 700.0
    },
    'iron': {
        'model': 'off',
        'solidus': 2200.0,
        'liquidus': 3000.0,
        'liquid_shear': 1.0e-5,
    }
}


find_partial_melter = ModelSearcher(viscosity, partial_melter_param_defaults)

# Helper Functions
@njit
def calc_melt_fraction(temperature: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the volumetric melt fraction

    :param temperature: <ndarray> Material/Layer temperature [K]
    :param solidus: <Float> Material/Layer solidus temperature [K]
    :param liquidus: <Float> Material/Layer liquidus temperature [K]
    :return: <Array> Volumetric Melt Fraction [m3 m-3]
    """

    melt_frac = (temperature - solidus) / (liquidus - solidus)

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    return melt_frac


@njit
def temp_from_melt(melt_frac: np.ndarray, solidus: float, liquidus: float) -> np.ndarray:
    """ Calculates the temperature from the volumetric melt fraction

    :param melt_frac: <ndarray> Volumetric Melt Fraction [m3 m-3]
    :param solidus: <Float> Material/Layer solidus temperature [K]
    :param liquidus: <Float> Material/Layer liquidus temperature [K]
    :return: <Array> Temperature [K]
    """

    # Check for over/under-shoots
    melt_frac[melt_frac < 0.] = 0.
    melt_frac[melt_frac > 1.] = 1.

    return melt_frac * (liquidus - solidus) + solidus



class PartialMelt(LayerModel):

    default_config = partial_melter_param_defaults
    config_key = 'partial_melt'


    def __init__(self, layer: LayerType):

        super().__init__(layer=layer, function_searcher=find_partial_melter, automate=True)

        self.calc_melt_fraction = calc_melt_fraction
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

        return self.func(self.layer.temperature, self.layer.melt_fraction, premelt_viscosity, premelt_shear, liquid_viscosity, *self.inputs)

    @staticmethod
    def _calculate_off(premelt_viscosity: np.ndarray, premelt_shear: np.ndarray, liquid_viscosity: np.ndarray):

        return premelt_viscosity, premelt_shear