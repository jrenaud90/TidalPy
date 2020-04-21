from typing import TYPE_CHECKING, Tuple

import numpy as np

from .complexCompliance import ComplexCompliance
from .defaults import rheology_defaults
from .partialMelt import PartialMelt
from .viscosity import SolidViscosity, LiquidViscosity
from .. import debug_mode
from ..exceptions import MissingArgumentError, ArgumentException, OuterscopeAttributeSetError, UnusualRealValueError
from ..initialize import log
from ..types import FloatArray
from ..utilities.arrayHelp import reshape_help
from ..utilities.classes import LayerConfigHolder

if TYPE_CHECKING:
    from ..structures.layers import ThermalLayer


class Rheology(LayerConfigHolder):

    """ Rheology class - Holder for all strength-related physics

    Rheology class stores model parameters and methods for viscosity(P, T, melt_frac), shear_modulus(P, T, melt_frac)
        and complex_compliance(viscosity, shear_modulus, forcing frequency).
    """

    default_config = rheology_defaults
    layer_config_key = 'rheology'

    def __init__(self, layer: 'ThermalLayer', store_config_in_layer: bool = True):

        super().__init__(layer, store_config_in_layer)

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

    def clear_state(self):

        super().clear_state()

        # Clear the state of all inner models
        for model in [self.viscosity_model, self.liquid_viscosity_model,
                      self.partial_melting_model, self.complex_compliance_model]:
            model.clear_state()


    def set_state(self, viscosity: FloatArray = None, shear_modulus: FloatArray = None, called_by_layer: bool = False):
        """ Set the rheology state and recalculate any parameters that may be affected.

        Setting the state manually overrides the functionality of the viscosity, shear, and partial melting functions.

        Parameters
        ----------
        viscosity : FloatArray
            Layer/Material viscosity [Pa s]
        shear_modulus : FloatArray
            Layer/Material shear modulus [Pa]
        called_by_layer : bool = False
            Should be False unless this function is called by the layer.set_strength method

        See Also
        --------
        TidalPy.structures.layers.layers.ThermalLayer.set_strength
        """

        # Now set the relevant parameters
        if viscosity is None and shear_modulus is None:
            ArgumentException('Function call has no effect.')
        if viscosity is None and self.viscosity is None:
            raise MissingArgumentError('Viscosity was not provided and is not already set.')
        if shear_modulus is None and self.shear_modulus is None:
            raise MissingArgumentError('Shear Modulus was not provided and is not already set.')

        # The shape of the viscosity and shear modulus must match the planet's global shape parameter (visco/shear.shape
        #    for each layer).
        if viscosity is not None:
            new_shape, viscosity = reshape_help(viscosity, self.world.global_shape,
                                                f'{self}.set_state.viscosity')
            if new_shape:
                self.world.change_shape(new_shape)
        if shear_modulus is not None:
            new_shape, shear_modulus = reshape_help(shear_modulus, self.world.global_shape,
                                                    f'{self}.set_state.shear_modulus')
            if new_shape:
                self.world.change_shape(new_shape)

        # Check for unusual values
        if debug_mode:
            for value in [viscosity, shear_modulus]:
                if value is None:
                    continue
                if np.any(value > 1.0e35) or np.any(value < 1.0e-10):
                    raise UnusualRealValueError

        # The premelt viscosity are no longer applicable as these new values would override them. To ensure they are
        #    not used let's set them to null
        self.liquid_viscosity_model._viscosity = None
        self.viscosity_model._viscosity = None

        # Now override the state variables
        if viscosity is not None:
            self.partial_melting_model._postmelt_viscosity = viscosity
        if shear_modulus is not None:
            self.partial_melting_model._postmelt_shear_modulus = shear_modulus
            self.partial_melting_model._postmelt_compliance = 1. / shear_modulus

        # A change to the viscosity or shear modulus will change the complex compliance so anything that depends on
        #     that will need to be updated as well.
        if not called_by_layer:
            self.world.update_tides()

    def update_strength(self, called_from_layer: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """ Calculates the strength of a layer/material based on temperature and pressure.

        Strength, in this context, refers to the layer's viscosity and shear modulus.

        Depending upon the model, partial melt fraction is also calculated and used to further modify the final
            strength.

        Parameters
        ----------
        called_from_layer: bool = False
            Should be false unless this method is being called from its host layer.

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

        # Now that there is a new viscosity and shear modulus, update complex compliance as well.
        #    Complex compliances depend upon forcing frequencies, which may or may not be set when this is called.
        if self.layer.is_tidal:
            if self.unique_tidal_frequencies is not None:
                self.complex_compliance_model.calculate()

        # A change to the viscosity or shear modulus will change the complex compliance so anything that depends on
        #     that will need to be updated as well.
        if not called_from_layer and self.layer.is_tidal:
            self.world.update_tides()

        return self.viscosity, self.shear_modulus


    # Inner-scope reference properties
    # # Viscosity Class
    @property
    def premelt_viscosity(self):
        return self.viscosity_model.viscosity

    @premelt_viscosity.setter
    def premelt_viscosity(self, value):
        self.viscosity_model.viscosity = value

    # # Liquid Viscosity Class
    @property
    def liquid_viscosity(self):
        return self.liquid_viscosity_model.viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        self.liquid_viscosity.liquid_viscosity = value

    # # Partial Melting Class
    @property
    def melt_fraction(self):
        return self.partial_melting_model.melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        self.partial_melting_model.melt_fraction = value

    @property
    def postmelt_viscosity(self):
        return self.partial_melting_model.postmelt_viscosity

    @postmelt_viscosity.setter
    def postmelt_viscosity(self, value):
        self.partial_melting_model.postmelt_viscosity = value

    @property
    def postmelt_shear_modulus(self):
        return self.partial_melting_model.postmelt_shear_modulus

    @postmelt_shear_modulus.setter
    def postmelt_shear_modulus(self, value):
        self.partial_melting_model.postmelt_shear_modulus = value

    @property
    def postmelt_compliance(self):
        return self.partial_melting_model.postmelt_compliance

    @postmelt_compliance.setter
    def postmelt_compliance(self, value):
        self.partial_melting_model.postmelt_compliance = value

    # # Complex Compliance Class
    @property
    def complex_compliances(self):
        return self.complex_compliance_model.complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        self.complex_compliance_model.complex_compliances = value


    # Outer-scope reference properties
    # # Layer Class
    @property
    def premelt_shear(self):
        return self.layer.static_shear_modulus

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise OuterscopeAttributeSetError

    @property
    def quality_factor(self):
        return self.layer.quality_factor

    @quality_factor.setter
    def quality_factor(self, value):
        raise OuterscopeAttributeSetError

    @property
    def beta(self):
        return self.layer.beta

    @beta.setter
    def beta(self, value):
        raise OuterscopeAttributeSetError

    # # Tides Class
    @property
    def tides(self):
        return self.world.tides

    @tides.setter
    def tides(self, value):
        raise OuterscopeAttributeSetError

    @property
    def unique_tidal_frequencies(self):
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