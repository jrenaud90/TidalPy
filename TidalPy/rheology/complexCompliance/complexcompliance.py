from typing import TYPE_CHECKING, List, Dict

import numpy as np

from ...exceptions import ImproperAttributeHandling, MissingAttributeError, AttributeNotSetError, IncorrectAttributeType
from ...utilities.model import LayerModelHolder
from . import known_models, known_model_live_args, known_model_const_args
from .defaults import complex_compliance_defaults


class ComplexCompliance(LayerModelHolder):

    """ Complex Compliance Class - Child of LayerModelHolder Class

    Complex compliance provides the functionality to calculate the complex compliance provided a viscosity, shear
        modulus, and forcing frequency. Different rheological models provide different functional forms of the
        complex compliance function.
    """

    default_config = complex_compliance_defaults
    known_models = known_models
    known_model_const_args = known_model_const_args
    known_model_live_args = known_model_live_args
    model_config_key = ('rheology', 'complex_compliance')

    def __init__(self, layer, rheology_class, model_name: str = None,
                 store_config_in_layer: bool = True):

        super().__init__(layer, model_name, store_config_in_layer)

        self.rheology_class = rheology_class
        self._complex_compliances = None

    def _calculate(self) -> Dict[str, np.ndarray]:
        """ Calculate the complex compliance for forcing frequency mode the planet is experiencing.

        Returns
        -------
        complex_compliance : List[np.ndarray]
            Complex compliance of the layer/material [Pa-1]
        """

        complex_compliances = {mode_name: self.func(tidal_freq, *self.live_inputs, *self.inputs)
                               for mode_name, tidal_freq in zip(self.tidal_mode_names, self.tidal_freqs)}

        self._complex_compliances = complex_compliances

        return complex_compliances

    # State properties
    @property
    def complex_compliances(self) -> Dict[str, np.ndarray]:
        return self._complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        raise ImproperAttributeHandling

    # Outerscope properties
    @property
    def tidal_freqs(self) -> List[np.ndarray]:
        return self.rheology_class.unique_tidal_freqs

    @tidal_freqs.setter
    def tidal_freqs(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_mode_names(self) -> List[str]:
        return self.rheology_class.tidal_mode_names

    @tidal_mode_names.setter
    def tidal_mode_names(self, value):
        raise ImproperAttributeHandling

    @property
    def compliance(self):
        return self.rheology_class.compliance

    @compliance.setter
    def compliance(self, value):
        raise ImproperAttributeHandling

    @property
    def viscosity(self):
        return self.rheology_class.viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise ImproperAttributeHandling

    @property
    def quality_factor(self):
        return self.layer.world.quality_factor

    @quality_factor.setter
    def quality_factor(self, value):
        raise ImproperAttributeHandling

    @property
    def beta(self):
        return self.layer.world.beta

    @beta.setter
    def beta(self, value):
        raise ImproperAttributeHandling