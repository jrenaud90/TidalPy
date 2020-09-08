from typing import Union

from .basic import BaseWorld
from .gas import GasGiantWorld
from .stellar import StarWorld
from .tidal import TidalWorld
from .layered import LayeredWorld
from .burnman import BurnManWorld

GasStarWorldType = Union[GasGiantWorld, StarWorld]
BurnManType = BurnManWorld
TidalWorldType = Union[GasStarWorldType, BurnManType, TidalWorld, LayeredWorld]

# FIXME: Add all world types here:
world_types = {
    'test': 1
    }