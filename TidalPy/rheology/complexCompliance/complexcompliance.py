from typing import TYPE_CHECKING

import numpy as np

from ...exceptions import ImproperAttributeHandling
from ...utilities.model import LayerModelHolder
from . import known_models, known_model_live_args, known_model_const_args
from .defaults import complex_compliance_defaults

if TYPE_CHECKING:
    from ...structures.layers import ThermalLayer
    from ..rheology import Rheology


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
    model_config_key = 'complex_compliance'

    def __init__(self, layer: ThermalLayer, rheology_class: Rheology, model_name: str = None,
                 store_config_in_layer: bool = True, call_reinit: bool = True):

        super().__init__(layer, model_name, store_config_in_layer, call_reinit)

        self.rheology_class = rheology_class

    def _calculate(self, compliance: np.ndarray, viscosity: np.ndarray, frequency: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        compliance : np.ndarray
            Layer/Material compliance (inverse of rigidity) [Pa-1]
        viscosity : np.ndarray
            Layer/Material effective (post-melting) viscosity [Pa s]
        frequency : np.ndarray
            Forcing frequency [rads s-1]
        other_inputs : tuple
            Constant and live arguments (needed for some functions)

        Returns
        -------
        complex_compliance : np.ndarray
            Complex compliance of the layer/material [Pa-1]
        """

        return self.func(compliance, viscosity, frequency, *self.inputs, *self.live_inputs)

    @property
    def zeta(self):
        return self.rheology_class.zeta

    @zeta.setter
    def zeta(self, value):
        raise ImproperAttributeHandling