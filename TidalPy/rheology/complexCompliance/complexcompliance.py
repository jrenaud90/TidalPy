from typing import List, Dict

import numpy as np

from . import known_models, known_model_live_args, known_model_const_args
from .defaults import complex_compliance_defaults
from ...exceptions import ImproperPropertyHandling
from ...tides.mode_manipulation import FreqSig
from ...utilities.classes.model import LayerModelHolder


class ComplexCompliance(LayerModelHolder):

    """ Complex Compliance Class - Child of LayerModelHolder Class

    Complex compliance provides the functionality to calculate the complex compliance once provided a viscosity, shear
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

    def clear_state(self):
        """ Clears the current state of the class without destroying data set by initialization.
        """

        super().clear_state()

        self._complex_compliances = None

    def _calculate(self) -> Dict[FreqSig, np.ndarray]:
        """ Calculate the complex compliance for forcing frequency mode the planet is experiencing.

        Returns
        -------
        complex_compliance : List[np.ndarray]
            Complex compliance of the layer/material [Pa-1]
        """

        complex_compliances = {mode_signature: self.func_array(tidal_freq, *self.live_inputs, *self.inputs)
                               for mode_signature, tidal_freq in self.tidal_freqs.items()}

        self._complex_compliances = complex_compliances

        return complex_compliances

    # State properties
    @property
    def complex_compliances(self) -> Dict[FreqSig, np.ndarray]:
        return self._complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        raise ImproperPropertyHandling

    # Outerscope properties
    @property
    def tidal_freqs(self):
        return self.rheology_class.unique_tidal_freqs

    @tidal_freqs.setter
    def tidal_freqs(self, value):
        raise ImproperPropertyHandling

    @property
    def compliance(self):
        return self.rheology_class.compliance

    @compliance.setter
    def compliance(self, value):
        raise ImproperPropertyHandling

    @property
    def viscosity(self):
        return self.rheology_class.viscosity

    @viscosity.setter
    def viscosity(self, value):
        raise ImproperPropertyHandling

    @property
    def quality_factor(self):
        return self.layer.world.quality_factor

    @quality_factor.setter
    def quality_factor(self, value):
        raise ImproperPropertyHandling

    @property
    def beta(self):
        return self.layer.world.beta

    @beta.setter
    def beta(self, value):
        raise ImproperPropertyHandling