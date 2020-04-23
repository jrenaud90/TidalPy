from typing import Union

from .basic import LayerBase
from .gas import GasLayer
from .thermal import ThermalLayer

LayerTypes = Union[ThermalLayer, LayerBase, GasLayer]