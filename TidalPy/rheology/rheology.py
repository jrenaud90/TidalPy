from typing import TYPE_CHECKING

import numpy as np

from .complex_compliance import ComplexCompliance
from .defaults import rheology_defaults
from .partial_melt import PartialMelt
from .viscosity import SolidViscosity, LiquidViscosity
from .. import debug_mode, log
from ..exceptions import MissingArgumentError, OuterscopePropertySetError, UnusualRealValueError, \
    InitiatedPropertyChangeError, IncorrectMethodToSetStateProperty
from ..utilities.classes import LayerConfigHolder
from ..utilities.types import FloatArray

if TYPE_CHECKING:
    from ..structures.layers import PhysicsLayer


class Rheology(LayerConfigHolder):

    """ Rheology
    A class for methods and attributes used to calculate a layer's strength-related physics.

    Rheology class stores model parameters and methods for viscosity(P, T, melt_frac), shear_modulus(P, T, melt_frac)
        and complex_compliance(viscosity, shear_modulus, forcing frequency).

    See Also
    --------
    TidalPy.utilities.methods.LayerConfigHolder
    TidalPy.rheology.complex_compliance.ComplexCompliance
    TidalPy.rheology.partial_melt.PartialMelt
    TidalPy.rheology.viscosity.ViscosityParentClass
    """

    default_config = rheology_defaults
    layer_config_key = 'rheology'

    def __init__(self, layer: 'PhysicsLayer', store_config_in_layer: bool = True):
        """ Constructor for Rheology class

        Parameters
        ----------
        layer : PhysicalLayerType
            The layer instance which the model should perform calculations on.
        store_config_in_layer: bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `layer.config` dictionary.
        """

        self.default_config = rheology_defaults[layer.type]

        super().__init__(layer, store_config_in_layer)

        # Initialized properties
        # Load in sub-models
        self._viscosity_model = \
            SolidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self._liquid_viscosity_model = \
            LiquidViscosity(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self._partial_melting_model = \
            PartialMelt(self.layer, self, store_config_in_layer=self.store_config_in_layer)
        self._complex_compliance_model = \
            ComplexCompliance(self.layer, self, store_config_in_layer=self.store_config_in_layer)

        # Report information about the models loaded
        log.debug(f"Rheology loaded into {self.layer}:\n"
                  f"\tSolid Viscosity:     {self.viscosity_model.model}\n"
                  f"\tLiquid Viscosity:    {self.liquid_viscosity_model.model}\n"
                  f"\tPartial Melting:     {self.partial_melting_model.model}\n"
                  f"\tComplex Compliance:  {self.complex_compliance_model.model}")

    def tidal_frequencies_changed(self, collapse_tidal_modes: bool = True):
        """ The tidal frequencies have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        log.debug(f'Tidal frequencies changed called for {self}.')

        if self.layer.is_tidal:
            if self.unique_tidal_frequencies is not None:
                # Calculate new complex compliances
                self.complex_compliance_model.calculate()

                # Tell the rheology class that the complex compliances have changed.
                self.complex_compliances_changed(collapse_tidal_modes=collapse_tidal_modes)

    def complex_compliances_changed(self, collapse_tidal_modes: bool = True):
        """ The complex compliances have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        log.debug(f'Complex compliances changed called for {self}.')

        if self.complex_compliances is not None:
            self.layer.complex_compliances_changed(collapse_tidal_modes=collapse_tidal_modes)

    def temperature_pressure_changed(self):
        """ The layer's temperature and/or pressure has changed. Make any necessary updates. """

        log.debug(f'Temperature change called for {self}.')

        if self.layer.temperature is not None:

            strength_changed = False
            # Update pre-partial melt solid viscosity
            if self.viscosity_model is not None:
                self.viscosity_model.calculate()
                strength_changed = True

            # Update liquid viscosity
            if self.liquid_viscosity_model is not None:
                self.liquid_viscosity_model.calculate()
                strength_changed = True

            # Update partial melting
            if self.partial_melting_model is not None:
                self.partial_melting_model.calculate()
                strength_changed = True

            if strength_changed:
                # Tell rheology that the strength has changed.
                self.strength_changed()

    def strength_changed(self):
        """ The layer's viscosity and/or shear modulus has changed. Make any necessary updates. """

        log.debug(f'Strength changed called for {self}.')

        if self.viscosity is not None and self.shear_modulus is not None:
            if self.complex_compliance_model is not None:
                # Calculate new complex compliances
                self.complex_compliance_model.calculate()

                # Tell the rheology class that the complex compliances have changed.
                self.complex_compliances_changed()

    def clear_state(self):
        """ Clear the state of all rheological parameters (all sub models' `clear_state` method will be called) """

        super().clear_state()

        # Clear the state of all inner models
        for model in [self.viscosity_model, self.liquid_viscosity_model,
                      self.partial_melting_model, self.complex_compliance_model]:
            model.clear_state()

    def set_state(self, viscosity: FloatArray = None, shear_modulus: FloatArray = None):
        """ Set the rheology state and recalculate any parameters that may be affected.

        Setting the state manually overrides the functionality of the viscosity, shear, and partial melting functions.

        Parameters
        ----------
        viscosity : FloatArray
            Layer/Material viscosity [Pa s]
        shear_modulus : FloatArray
            Layer/Material shear modulus [Pa]

        See Also
        --------
        TidalPy.structures.layers.PhysicsLayer.set_strength
        """

        strength_changed = False

        # Now set the relevant parameters
        if viscosity is None and shear_modulus is None:
            log.debug('Rheology set_state called but no arguments were provided.')

        # It is fine to only pass one state property, but we need the other to perform calculations. The method will
        #    check if the other property is already set, but if it isn't then an exception will be raised.
        if viscosity is None and self.viscosity is None:
            log.error('Rheology set_state called and viscosity is not set and was not provided.')
            raise MissingArgumentError('Viscosity was not provided and is not already set.')
        if shear_modulus is None and self.shear_modulus is None:
            log.error('Rheology set_state called and shear modulus is not set and was not provided.')
            raise MissingArgumentError('Shear Modulus was not provided and is not already set.')

        # Check for unusual values
        if debug_mode:
            for value in [viscosity, shear_modulus]:
                if value is None:
                    continue
                if np.any(value > 1.0e35) or np.any(value < 1.0e-10):
                    raise UnusualRealValueError

        # The premelt strength is no longer applicable as these new values would override them. To ensure they are
        #    not used let's set them to null
        if viscosity is not None:
            self.liquid_viscosity_model._viscosity = None
            self.viscosity_model._viscosity = None

        # Now override the state properties
        if viscosity is not None:
            self.partial_melting_model._postmelt_viscosity = viscosity
            strength_changed = True
        if shear_modulus is not None:
            self.partial_melting_model._postmelt_shear_modulus = shear_modulus
            self.partial_melting_model._postmelt_compliance = 1. / shear_modulus
            strength_changed = True

        if strength_changed:
            # Tell rheology that the strength has changed.
            self.strength_changed()


    # # Initialized properties
    @property
    def viscosity_model(self) -> SolidViscosity:
        """ Viscosity model instance used to calculate layer's solid viscosity (no partial melt) """
        return self._viscosity_model

    @viscosity_model.setter
    def viscosity_model(self, value):
        raise InitiatedPropertyChangeError

    @property
    def liquid_viscosity_model(self) -> LiquidViscosity:
        """ Viscosity model instance used to calculate layer's liquid viscosity (no partial melt) """
        return self._liquid_viscosity_model

    @liquid_viscosity_model.setter
    def liquid_viscosity_model(self, value):
        raise InitiatedPropertyChangeError

    @property
    def partial_melting_model(self) -> PartialMelt:
        """ Partial melting model instance used to modify the solid & liquid viscosities based on the melt fraction """
        return self._partial_melting_model

    @partial_melting_model.setter
    def partial_melting_model(self, value):
        raise InitiatedPropertyChangeError

    @property
    def complex_compliance_model(self) -> ComplexCompliance:
        """ Complex compliance model instance used to calculate the layers complex compliance """
        return self._complex_compliance_model

    @complex_compliance_model.setter
    def complex_compliance_model(self, value):
        raise InitiatedPropertyChangeError


    # Inner-scope reference properties
    # # Viscosity Class
    @property
    def premelt_viscosity(self):
        """ Inner-scope wrapper for sold_viscosity_model.viscosity """
        return self.viscosity_model.viscosity

    @premelt_viscosity.setter
    def premelt_viscosity(self, value):
        raise IncorrectMethodToSetStateProperty

    # # Liquid Viscosity Class
    @property
    def liquid_viscosity(self):
        """ Inner-scope wrapper for liquid_viscosity_model.viscosity """
        return self.liquid_viscosity_model.viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        raise IncorrectMethodToSetStateProperty

    # # Partial Melting Class
    @property
    def melt_fraction(self):
        """ Inner-scope wrapper for partial_melting_model.melt_fraction """
        return self.partial_melting_model.melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def postmelt_viscosity(self):
        """ Inner-scope wrapper for partial_melting_model.postmelt_viscosity """
        return self.partial_melting_model.postmelt_viscosity

    @postmelt_viscosity.setter
    def postmelt_viscosity(self, postmelt_viscosity: FloatArray):
        self.set_state(viscosity=postmelt_viscosity)

    @property
    def postmelt_shear_modulus(self):
        """ Inner-scope wrapper for partial_melting_model.postmelt_shear_modulus """
        return self.partial_melting_model.postmelt_shear_modulus

    @postmelt_shear_modulus.setter
    def postmelt_shear_modulus(self, postmelt_shear_modulus: FloatArray):
        self.set_state(shear_modulus=postmelt_shear_modulus)

    @property
    def postmelt_compliance(self):
        """ Inner-scope wrapper for partial_melting_model.postmelt_compliance """
        return self.partial_melting_model.postmelt_compliance

    @postmelt_compliance.setter
    def postmelt_compliance(self, postmelt_compliance: FloatArray):
        self.set_state(shear_modulus=postmelt_compliance**(-1))

    # # Complex Compliance Class
    @property
    def complex_compliances(self):
        """ Inner-scope wrapper for complex_compliance_model.complex_compliances """
        return self.complex_compliance_model.complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        raise IncorrectMethodToSetStateProperty


    # Outer-scope reference properties
    # # Layer Class
    @property
    def premelt_shear(self):
        """ Outer-scope wrapper for layer.static_shear_modulus """
        return self.layer.static_shear_modulus

    @premelt_shear.setter
    def premelt_shear(self, value):
        raise OuterscopePropertySetError

    @property
    def quality_factor(self):
        """ Outer-scope wrapper for layer.quality_factor """
        return self.layer.quality_factor

    @quality_factor.setter
    def quality_factor(self, value):
        raise OuterscopePropertySetError

    @property
    def beta(self):
        """ Outer-scope wrapper for layer.beta """
        return self.layer.beta

    @beta.setter
    def beta(self, value):
        raise OuterscopePropertySetError

    # World Class
    @property
    def unique_tidal_frequencies(self):
        """ Outer-scope wrapper for world.unique_tidal_frequencies """
        return self.world.unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise OuterscopePropertySetError


    # # Aliased properties
    @property
    def viscosity(self):
        """ Alias for self.postmelt_viscosity """
        return self.postmelt_viscosity

    @viscosity.setter
    def viscosity(self, value):
        self.postmelt_viscosity = value

    @property
    def shear_modulus(self):
        """ Alias for self.postmelt_shear_modulus """
        return self.postmelt_shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        self.postmelt_shear_modulus = value

    @property
    def shear(self):
        """ Alias for self.postmelt_shear_modulus """
        return self.postmelt_shear_modulus

    @shear.setter
    def shear(self, value):
        self.postmelt_shear_modulus = value

    @property
    def compliance(self):
        """ Alias for self.postmelt_compliance """
        return self.postmelt_compliance

    @compliance.setter
    def compliance(self, value):
        self.postmelt_compliance = value
