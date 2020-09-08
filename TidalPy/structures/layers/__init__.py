from typing import Union

from .basic import LayerBase
from .gas import GasLayer
from .physics import PhysicsLayer

LayerTypes = Union[PhysicsLayer, LayerBase, GasLayer]