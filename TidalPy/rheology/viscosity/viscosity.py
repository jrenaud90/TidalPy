from typing import TYPE_CHECKING

import numpy as np

from TidalPy.exceptions import (IncorrectMethodToSetStateProperty, InitiatedPropertyChangeError,
                                OuterscopePropertySetError)
from TidalPy.utilities.classes.model import LayerModelHolder

from . import known_model_const_args, known_model_live_args, known_models
if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray
    from TidalPy.rheology import Rheology
    from TidalPy.structures.layers import PhysicalLayerType


class ViscosityParentClass(LayerModelHolder):
    """ ViscosityParentClass

    Common Viscosity class for liquid and solid viscosity models.

    See Also
    --------
    TidalPy.utilities.methods.model.LayerModelHolder
    TidalPy.rheology.Rheology
    TidalPy.rheology.viscosity.LiquidViscosityParent
    TidalPy.rheology.viscosity.SolidViscosityParent
    """

    is_liquid = False

    def __init__(
        self, layer: 'PhysicalLayerType', rheology_class: 'Rheology', model_name: str = None,
        store_config_in_layer: bool = True, initialize: bool = True
        ):
        """ Constructor for Parent Viscosity class

        Common Viscosity class for liquid and solid viscosity models.

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

        # State properties
        self._viscosity = None

        if initialize:
            self.reinit(initial_init=True)

    def clear_state(self):

        super().clear_state()

        self._viscosity = None

    def _calculate(self) -> 'FloatArray':
        """ Wrapper for the viscosity calculator

        Returns
        -------
        viscosity : FloatArray
            Viscosity of the layer [Pa s]
        """

        if type(self.temperature) == np.ndarray:
            viscosity = self.func_array(*self.live_inputs, *self.inputs)
        else:
            viscosity = self.func(*self.live_inputs, *self.inputs)

        self._viscosity = viscosity

        return viscosity

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
    def viscosity(self) -> 'FloatArray':
        """ Viscosity (Solid or Liquid depending on the class type) no partial melting effects have been applied """
        return self._viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise IncorrectMethodToSetStateProperty

    # # Outer-scope references
    @property
    def temperature(self):
        """ Outer-scope wrapper for layer.temperature """
        return self.layer.temperature

    @temperature.setter
    def temperature(self, value):
        raise OuterscopePropertySetError

    @property
    def pressure(self) -> 'FloatArray':
        """ Outer-scope wrapper for layer.pressure """
        # Pressure is often turned off in various model runs, so it may not be set. If that is the case let's make sure
        #    that zeros are provided.
        _pressure = self.layer.pressure
        if _pressure is None or not self.layer.use_pressure_in_strength_calc:
            _pressure = np.zeros_like(self.temperature)

        return _pressure

    @pressure.setter
    def pressure(self, value):
        raise OuterscopePropertySetError


class LiquidViscosity(ViscosityParentClass):
    """ LiquidViscosity
    Model for calculating viscosity of a liquid.

    See Also
    --------
    TidalPy.utilities.methods.model.LayerModelHolder
    TidalPy.rheology.Rheology
    TidalPy.rheology.viscosity.ViscosityParentClass
    """

    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = 'liquid_viscosity'
    is_liquid = True


class SolidViscosity(ViscosityParentClass):
    """ SolidViscosity
    Model for calculating viscosity of a solid.


    See Also
    --------
    TidalPy.utilities.methods.model.LayerModelHolder
    TidalPy.rheology.Rheology
    TidalPy.rheology.viscosity.ViscosityParentClass
    """

    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = 'solid_viscosity'
    is_liquid = False
