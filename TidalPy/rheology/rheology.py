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
from ..tides.love1d import effective_rigidity_general, complex_love_general


class Rheology(LayerConfigHolder):

    """ Rheology class - Holder for all strength-related physics

    Rheology class stores model parameters and methods for viscosity(P, T, melt_frac), shear_modulus(P, T, melt_frac)
        and complex_compliance(viscosity, shear_modulus, forcing frequency).
    """

    default_config = rheology_defaults
    layer_config_key = 'rheology'

    def __init__(self, layer,
                 orbital_truncation_level: int = 2, tidal_order_lvl: int = 2, use_nsr: bool = True,
                 store_config_in_layer: bool = True):

        super().__init__(layer, store_config_in_layer)

        # State properties
        self._effective_rigidities = None
        self._complex_love_numbers = None

        # Pull out model information
        self.order_lvl = tidal_order_lvl

        # Pull out switches
        self.use_planet_params_for_love_calc = self.config['use_planet_params_for_love_calc']

        # Pull out parameters for Love number and effective rigidity calculations
        if self.use_planet_params_for_love_calc:
            self.radius = self.layer.world.radius
            self.density = self.layer.world.density
            self.gravity = self.layer.world.gravity
        else:
            self.radius = self.layer.radius
            self.density = self.layer.density_bulkd
            self.gravity = self.layer.gravity_surf

        self.order_l_list = range(2, self.order_lvl+1)

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
            f"\tSolid Viscosity:    {self.viscosity_model.model}\n"
            f"\tLiquid Viscosity:   {self.liquid_viscosity_model.model}\n"
            f"\tPartial Melting:    {self.partial_melting_model.model}\n"
            f"\tComplex Compliance: {self.complex_compliance_model.model}", level='info')

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

    def calculate_love(self, only_frequency_change: bool = False) -> \
            Tuple[Dict[str, np.ndarray], List[Dict[str, np.ndarray]]]:
        """ Calculate the complex compliances and Love numbers based on the post-melting viscosity and shear modulus

        Returns
        -------
        complex_compliances : Dict[str, np.ndarray]
            Complex compliances in a dictionary for each tidal frequency [Pa-1]
        complex_love_numbers : List[Dict[str, np.ndarray]
            Complex love numbers stored in a dictionary for each tidal frequency
            Multiple dictionaries are stored in a list over tidal order l's
        """

        # Calculate complex compliances for each tidal frequency
        self.complex_compliance_model.calculate()

        if not only_frequency_change:
            # Calculate effective rigidity at each order l (effective rigidity does not change with frequency changes).
            self._effective_rigidities = [effective_rigidity_general(self.shear, self.gravity, self.radius, self.density, l)
                                          for l in self.order_l_list]

        # Calculate the complex love number at each tidal frequency and each order l
        self._complex_love_numbers = [{mode: complex_love_general(cmplx_comp, self.shear, eff_rigid, l)
                                       for mode, cmplx_comp in self.complex_compliances.items()}
                                      for eff_rigid, l in zip(self.effective_rigidities, self.order_l_list)]

        return self.complex_compliances, self.complex_love_numbers

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
    @property
    def premelt_shear(self):
        return self.layer.static_shear_modulus

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise OuterscopeAttributeSetError

    @property
    def tidal_modes(self):
        return self.layer.tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value):
        raise OuterscopeAttributeSetError

    @property
    def unique_tidal_freqs(self) -> List[np.ndarray]:
        return self.layer.unique_tidal_freqs

    @unique_tidal_freqs.setter
    def unique_tidal_freqs(self, value):
        raise OuterscopeAttributeSetError

    @property
    def tidal_mode_names(self) -> List[str]:
        return self.layer.tidal_mode_names

    @tidal_mode_names.setter
    def tidal_mode_names(self, value):
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