from typing import Union

from .basic import BaseWorld
from .burnman import BurnManWorld
from .gas import GasGiantWorld
from .layered import LayeredWorld
from .stellar import StarWorld
from .tidal import TidalWorld

GasStarWorldType = Union[GasGiantWorld, StarWorld]
BurnManType = BurnManWorld
TidalWorldType = Union[GasStarWorldType, BurnManType, TidalWorld, LayeredWorld]

# FIXME: Add all world types here:
world_types = {
    'star': StarWorld,
    'gas': GasGiantWorld,
    'gas_giant': GasGiantWorld,
    'ice_giant': GasGiantWorld,
    'simple_tide': TidalWorld,
    'layered': LayeredWorld,
    'burnman': BurnManWorld
}