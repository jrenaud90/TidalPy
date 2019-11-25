from typing import TYPE_CHECKING

import numpy as np

from ...exceptions import BadValueError
from ...performance import njit
from ...utilities.model import LayerModelHolder
from . import known_model_live_args, known_model_const_args, known_models
from .defaults import partial_melt_defaults

if TYPE_CHECKING:
    from ...structures.layers import ThermalLayer
    from ..rheology import Rheology


@njit
def calc_partial_melting(temperature: np.ndarray, solidus: float, liquidus: float):
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
    model_config_key = 'partial_melting'

    def __init__(self, layer: ThermalLayer, rheology_class: Rheology, model_name: str = None,
                 store_config_in_layer: bool = True, call_reinit: bool = True):

        super().__init__(layer, model_name, store_config_in_layer, call_reinit)

        self.rheology_class = rheology_class

        # Pull out specific information related to this module
        self.solidus = self.get_param('solidus')
        self.liquidus = self.get_param('liquidus')
        self.use_partial_melt = True
        if self.model == 'off':
            self.use_partial_melt = False

        # TODO: No debug calculation has been implemented
        self._calculate_debug = self._calculate


    def _calculate(self, temperature: np.ndarray, premelt_viscosity: np.ndarray, premelt_shear: np.ndarray,
                   liquid_viscosity: np.ndarray, other_inputs: tuple = tuple()):
        """ Wrapper for the partial melting function.

        First the partial melt will be updated based on the provided temperature. Then changes to the shear modulus
            and viscosity are updated.

        Parameters
        ----------
        temperature : np.ndarray
            Layer/Material temperature [K]
        premelt_viscosity : np.ndarray
            Layer/Material viscosity before partial melting is considered [Pa s]
        premelt_shear : np.ndarray
            Layer/Material shear modulus before partial melting is considered [Pa]
        liquid_viscosity : np.ndarray
            Layer/Material viscosity if it were completely molten at this temperature [Pa s]
        other_inputs : tuple
            Constant and live arguments (needed for some functions)

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
            melt_fraction = calc_partial_melting(temperature, self.solidus, self.liquidus)
            postmelt_viscosity, postmelt_shear_modulus = self.func(self.layer.temperature, self.layer.melt_fraction,
                                                                   premelt_viscosity, premelt_shear, liquid_viscosity,
                                                                   *self.inputs, *self.live_inputs)
        else:
            melt_fraction = np.zeros_like(temperature)
            postmelt_viscosity, postmelt_shear_modulus = premelt_viscosity, premelt_shear

        return melt_fraction, postmelt_viscosity, postmelt_shear_modulus