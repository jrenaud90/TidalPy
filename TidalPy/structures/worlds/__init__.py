from typing import Union

from .basic import GeometricWorld, WorldBase
from .gas import GasGiant
from .stellar import Star
from .tidal import TidalWorld, SimpleTidalWorld, LayeredWorld

TidalWorldType = Union[SimpleTidalWorld, TidalWorld]

# FIXME: Add all world types here:
world_types = {
    'test': 1
    }