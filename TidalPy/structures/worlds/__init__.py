from typing import Union

from .basic import GeometricWorld, WorldBase
from .gas import GasGiant
from .stellar import Star
from .tidal import TidalWorld, SimpleTidalWorld

TidalWorldType = Union[SimpleTidalWorld, TidalWorld]