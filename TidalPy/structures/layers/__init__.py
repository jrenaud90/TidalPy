from typing import Union

from .basic import LayerBase
from .gas import GasLayer
from .physics import PhysicsLayer

LayerType = Union[PhysicsLayer, LayerBase, GasLayer]
PhysicalLayerType = PhysicsLayer

known_layer_classes = {
    'gas'    : GasLayer,
    'physics': PhysicsLayer
    }

layers_class_by_world_class = {
    'layered'          : PhysicsLayer,
    'gas_giant_layered': GasLayer,
    }
