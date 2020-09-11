from typing import Union

from .basic import LayerBase
from .burnman import BurnManLayer
from .gas import GasLayer
from .physics import PhysicsLayer

LayerType = Union[PhysicsLayer, LayerBase, GasLayer, BurnManLayer]

known_layer_classes = {
    'gas': GasLayer,
    'physics': PhysicsLayer,
    'burnman': BurnManLayer
}

layers_class_by_world_class = {
    'burnman': BurnManLayer,
    'layered': PhysicsLayer,
    'gas_giant_layered': GasLayer,
}