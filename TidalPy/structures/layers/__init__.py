from typing import Union

from .thermal import ThermalLayer
from .basic import BasicLayer

LayerTypes = Union[ThermalLayer, BasicLayer, GasLayer]