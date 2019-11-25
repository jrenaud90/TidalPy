from typing import TYPE_CHECKING, List

import numpy as np

from ..types import ArrayNone
from ..utilities.classes import LayerConfigHolder
from .complexCompliance import ComplexCompliance
from .partialMelt import PartialMelt
from .viscosity import SolidViscosity, LiquidViscosity
from .defaults import rheology_defaults

if TYPE_CHECKING:
    from ..structures.layers import ThermalLayer


class Rheology(LayerConfigHolder):

    default_config = rheology_defaults
    layer_config_key = 'rheology'

    def __init__(self, layer: ThermalLayer, orbital_truncation_level: int = 2, tidal_order_l: int = 2,
                 store_config_in_layer: bool = True):

        super().__init__(layer, store_config_in_layer)

        # Pull out model information
        self.order_l = tidal_order_l
        self.orbit_trunc = orbital_truncation_level

        # State Attributes (live attributes required by some sub-modules)
        self._zeta = None

        # Load in sub-modules
        self.complex_compliance = ComplexCompliance(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.partial_melting = PartialMelt(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.viscosity = SolidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.liquid_viscosity = LiquidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.andrade_frequency = ...

    def calculate_strength(self, temperature: np.ndarray, pressure: ArrayNone = None, change_state: bool = True):
        """ Calculates the strength of a layer/material based on temperature and pressure.

        Strength, in this context, refers to the layer's viscosity and shear modulus.

        Depending upon the model, partial melt fraction is also calculated and used to further modify the final
            strength.

        Parameters
        ----------
        temperature : np.ndarray
            Layer temperature [K]
        pressure : ArrayNone = None
            Layer pressure [Pa]
            This can be set to None if ignoring pressure-dependence.
        change_state : bool = True

        Returns
        -------
        postmelt_viscosity : np.ndarray
            Final viscosity [Pa s]
        postmelt_shear_modulus : np.ndarray
            Final shear modulus [Pa]
        """

        # Calculate the pre-melting viscosities
        premelt_solid_viscosity = self.viscosity.calculate(temperature, pressure)
        liquid_viscosity = self.liquid_viscosity.calculate(temperature, pressure)
        premelt_shear_modulus = self.layer.static_shear_modulus

        # Add in the effects of partial melting
        melt_fraction, postmelt_viscosity, postmelt_shear_modulus = \
            self.partial_melting.calculate(temperature, premelt_solid_viscosity,
                                           premelt_shear_modulus, liquid_viscosity)

        if change_state and self.layer is not None:
            self.layer._melt_fraction = melt_fraction
            self.layer._viscosity = postmelt_viscosity
            self.layer._shear_modulus = postmelt_shear_modulus
            self.layer._compliance = postmelt_shear_modulus**(-1)

        return postmelt_viscosity, postmelt_shear_modulus


    def calculate_rheo






    # Live attributes required by some sub-modules
    @property
    def zeta(self) -> float:
        return self._zeta

    @zeta.setter
    def zeta(self, value: float):
        self._zeta = value