from typing import TYPE_CHECKING, List, Dict, Tuple
from functools import partial

import numpy as np

from ..exceptions import (ImproperAttributeHandling, ParameterValueError, MissingArgumentError, ArgumentOverloadError,
                          ArgumentException, ImplementationException, OuterscopeAttributeSetError)
from ..utilities.classes import LayerConfigHolder
from ..initialize import log
from .complexCompliance import ComplexCompliance
from .partialMelt import PartialMelt
from .viscosity import SolidViscosity, LiquidViscosity
from .defaults import rheology_defaults
from ..types import FloatArray

if TYPE_CHECKING:
    from ..structures.layers import ThermalLayer
    from ..tides.tides import Tides


class Rheology(LayerConfigHolder):

    """ Rheology class - Holder for all strength-related physics

    Rheology class stores model parameters and methods for viscosity(P, T, melt_frac), shear_modulus(P, T, melt_frac)
        and complex_compliance(viscosity, shear_modulus, forcing frequency).
    """

    default_config = rheology_defaults
    layer_config_key = 'rheology'

    def __init__(self, layer: 'ThermalLayer', store_config_in_layer: bool = True):

        super().__init__(layer, store_config_in_layer)

        # State properties
        self._effective_rigidities = None
        self._complex_love_numbers = None

        # Pull out switches
        self.use_planet_params_for_love_calc = self.config['use_planet_params_for_love_calc']

        # Load in sub-modules
        self.viscosity_model = \
            SolidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.liquid_viscosity_model = \
            LiquidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.partial_melting_model = \
            PartialMelt(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self.complex_compliance_model = \
            ComplexCompliance(self.layer, self, store_config_in_layer=self.store_config_in_layer)

        # Report information about the models loaded
        log(f"Rheology loaded into {self.layer}:\n"
            f"    Solid Viscosity:    {self.viscosity_model.model}\n"
            f"    Liquid Viscosity:   {self.liquid_viscosity_model.model}\n"
            f"    Partial Melting:    {self.partial_melting_model.model}\n"
            f"    Complex Compliance: {self.complex_compliance_model.model}", level='info')

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

        # A change to the viscosity or shear modulus will change the complex compliance so it needs to be updated
        self.calculate_love()

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

    def update_tides(self):
        """ Tell the tide class that changes have been made to the rheology that will impact tidal calculations.

        The tides module will update just the results based on this layer.

        """


        self.tides.rheology_update(self.layer)


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


    # State properties
    @property
    def effective_rigidities(self) -> List[np.ndarray]:
        return self._effective_rigidities

    @effective_rigidities.setter
    def effective_rigidities(self, value):
        raise ImproperAttributeHandling

    @property
    def complex_love_numbers(self) -> List[Dict[str, np.ndarray]]:
        return self._complex_love_numbers

    @complex_love_numbers.setter
    def complex_love_numbers(self, value):
        raise ImproperAttributeHandling


    # Outerscope reference properties
    #    Layer properties
    @property
    def premelt_shear(self):
        return self.layer.static_shear_modulus

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise OuterscopeAttributeSetError

    @property
    def quality_factor(self) -> float:
        return self.layer.quality_factor

    @quality_factor.setter
    def quality_factor(self, value):
        raise OuterscopeAttributeSetError

    @property
    def beta(self) -> float:
        return self.layer.beta

    @beta.setter
    def beta(self, value):
        raise OuterscopeAttributeSetError

    #    Tides properties
    @property
    def tides(self) -> Tides:
        return self.world.tides

    @tides.setter
    def tides(self, value):
        raise OuterscopeAttributeSetError

    @property
    def unique_tidal_frequencies(self) -> Dict[str, FloatArray]:
        return self.world.tides.unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise OuterscopeAttributeSetError


    # Alias properties
    @property
    def viscosity(self):
        return self.postmelt_viscosity

    @viscosity.setter
    def viscosity(self, value):
        self.postmelt_viscosity = value

    @property
    def shear_modulus(self):
        return self.postmelt_shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        self.postmelt_shear_modulus = value

    @property
    def shear(self):
        return self.postmelt_shear_modulus

    @shear.setter
    def shear(self, value):
        self.postmelt_shear_modulus = value

    @property
    def compliance(self):
        return self.postmelt_compliance

    @compliance.setter
    def compliance(self, value):
        self.postmelt_compliance = value