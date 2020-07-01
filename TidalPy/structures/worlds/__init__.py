from typing import Union

from .basic import GeometricWorld, BaseWorld
from .gas import GasGiantWorld
from .stellar import StarWorld
from .tidal import TidalWorld, SimpleTidalWorld
from .layered import LayeredWorld

GasPlasmaWorldType = Union[GasGiantWorld, StarWorld]
TidalWorldType = Union[GasPlasmaWorldType, TidalWorld, LayeredWorld]

# FIXME: Add all world types here:
world_types = {
    'test': 1
    }