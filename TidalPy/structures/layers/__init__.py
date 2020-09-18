from typing import Union

from .basic import LayerBase
from .burnman import BurnmanLayer
from .gas import GasLayer
from .physics import PhysicsLayer

LayerType = Union[PhysicsLayer, LayerBase, GasLayer, BurnmanLayer]
PhysicalLayerType = Union[PhysicsLayer, BurnmanLayer]

known_layer_classes = {
    'gas': GasLayer,
    'physics': PhysicsLayer,
    'burnman': BurnmanLayer
}

layers_class_by_world_class = {
    'burnman': BurnmanLayer,
    'layered': PhysicsLayer,
    'gas_giant_layered': GasLayer,
}