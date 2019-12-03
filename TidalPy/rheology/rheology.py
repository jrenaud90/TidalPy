from typing import TYPE_CHECKING, List
from functools import partial

import numpy as np

from ..dynamics import mode_types
from ..exceptions import (ImproperAttributeHandling, ParameterValueError, MissingArgumentError, ArgumentOverloadError,
                          ArgumentException)
from ..types import ArrayNone
from ..utilities.classes import LayerConfigHolder
from ..tides import calculate_tides
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
                 use_nsr: bool = True,
                 store_config_in_layer: bool = True):

        super().__init__(layer, store_config_in_layer)

        # Pull out model information
        self.order_l = tidal_order_l
        self.orbit_trunc_lvl = orbital_truncation_level
        self.use_nsr = use_nsr

        # Find the correct Love and effective-rigidity functions based on the order number.
        if self.order_l > 3:
            raise NotImplementedError(f'Tidal order {self.order_l} has not been implemented yet.')
        if self.orbit_trunc_lvl % 2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self.orbit_trunc_lvl <= 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self.orbit_trunc_lvl not in [2, 4, 6]:
            raise NotImplementedError(f'Orbital truncation level of {self.orbit_trunc_lvl} is not currently supported.')

        # Load in sub-modules
        self.complex_compliance_model = \
            ComplexCompliance(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.partial_melting_model = \
            PartialMelt(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.viscosity_model = \
            SolidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.liquid_viscosity_model = \
            LiquidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)

    def set_state(self, viscosity: np.ndarray = None, shear_modulus: np.ndarray = None):
        """ Set the rheology state and recalculate any parameters that may now be different.

        Parameters
        ----------
        viscosity : np.ndarray
            Layer/Material viscosity [Pa s]
        shear_modulus : np.ndarray
            Layer/Material shear modulus [Pa]
        """

        # The premelt viscosity are no longer applicable as these new values would override them. To ensure they are
        #    not used lets reset them
        self.liquid_viscosity_model._viscosity = None
        self.viscosity_model._viscosity = None

        # Now set the relevant parameters
        if viscosity is None and shear_modulus is None:
            ArgumentException('Function call has no effect.')
        if viscosity is None and self.viscosity is None:
            raise MissingArgumentError('Viscosity was not provided and is not already set.')
        if shear_modulus is None and self.shear_modulus is None:
            raise MissingArgumentError('Shear Modulus was not provided and is not already set.')

        # Now override the state variables
        if viscosity is not None:
            self.partial_melting_model._postmelt_viscosity = viscosity
        if shear_modulus is not None:
            self.partial_melting_model._postmelt_shear_modulus = shear_modulus
            self.partial_melting_model._postmelt_compliance = shear_modulus**(-1)

    def calculate_strength(self):
        """ Calculates the strength of a layer/material based on temperature and pressure.

        Strength, in this context, refers to the layer's viscosity and shear modulus.

        Depending upon the model, partial melt fraction is also calculated and used to further modify the final
            strength.

        Returns
        -------
        postmelt_viscosity : np.ndarray
            Final viscosity [Pa s]
        postmelt_shear_modulus : np.ndarray
            Final shear modulus [Pa]
        """

        # Calculate the pre-melting viscosities
        self.viscosity_model.calculate()
        self.liquid_viscosity_model.calculate()

        # Add in the effects of partial melting
        self.partial_melting_model.calculate()

        return self.viscosity, self.shear_modulus

    def calculate_compliances(self) -> List[np.ndarray]:
        """ Calculate the complex compliances based on the post-melting viscosity and shear modulus


        Returns
        -------
        complex_compliances : List[np.ndarray]
        """

        self.complex_compliance_model.calculate()

        return self.complex_compliances

    # State properties

    # Alias properties
    @property
    def viscosity(self):
        return self.postmelt_viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def shear_modulus(self):
        return self.postmelt_shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        raise ImproperAttributeHandling

    @property
    def shear(self):
        return self.postmelt_shear_modulus

    @shear.setter
    def shear(self, value):
        raise ImproperAttributeHandling

    # Innerscope reference properties
    @property
    def premelt_viscosity(self):
        return self.viscosity_model.viscosity

    @premelt_viscosity.setter
    def premelt_viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def liquid_viscosity(self):
        return self.liquid_viscosity_model.viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def melt_fraction(self):
        return self.partial_melting_model.melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        raise ImproperAttributeHandling

    @property
    def postmelt_viscosity(self):
        return self.partial_melting_model.postmelt_viscosity

    @postmelt_viscosity.setter
    def postmelt_viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def postmelt_shear_modulus(self):
        return self.partial_melting_model.postmelt_shear_modulus

    @postmelt_shear_modulus.setter
    def postmelt_shear_modulus(self, value):
        raise ImproperAttributeHandling

    @property
    def postmelt_compliance(self):
        return self.partial_melting_model.postmelt_compliance

    @postmelt_compliance.setter
    def postmelt_compliance(self, value):
        raise ImproperAttributeHandling

    @property
    def complex_compliances(self):
        return self.complex_compliance_model.complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        raise ImproperAttributeHandling

    # Outerscope reference properties
    @property
    def premelt_shear(self):
        return self.layer.static_shear_modulus

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_modes(self):
        return self.layer.world.tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value):
        raise ImproperAttributeHandling

    @property
    def quality_factor(self) -> float:
        return self.layer.world.quality_factor

    @quality_factor.setter
    def quality_factor(self, value):
        raise ImproperAttributeHandling

    @property
    def beta(self) -> float:
        return self.layer.beta

    @beta.setter
    def beta(self, value):
        raise ImproperAttributeHandling