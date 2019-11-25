from typing import TYPE_CHECKING

from ..utilities.classes import LayerConfigHolder
from .complexCompliance import ComplexCompliance
from .partialMelt import PartialMelt
from .defaults import rheology_defaults

if TYPE_CHECKING:
    from ..structures.layers import ThermalLayer


class Rheology(LayerConfigHolder):

    default_config = rheology_defaults
    layer_config_key = 'rheology'

    def __init__(self, layer: ThermalLayer, store_config_in_layer: bool = True, call_reinit: bool = True):

        super().__init__(layer, store_config_in_layer, call_reinit)

        # Load in sub-modules
        self.complex_compliance = ComplexCompliance(self.layer, self, store_config_in_layer=store_config_in_layer,
                                                    call_reinit=call_reinit)
        self.partial_melting = PartialMelt(self.layer, self, store_config_in_layer=store_config_in_layer,
                                           call_reinit=call_reinit)
        self.viscosity =