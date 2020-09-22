from typing import Union

from .basic import BaseWorld
from .burnman import BurnManWorld
from .gas import GasGiantWorld, GasGiantLayeredWorld
from .layered import LayeredWorld
from .stellar import StarWorld
from .tidal import TidalWorld

AllWorldType = Union[BaseWorld, TidalWorld, GasGiantWorld, StarWorld, LayeredWorld, GasGiantLayeredWorld, BurnManWorld]
all_world_types = (BaseWorld, TidalWorld, GasGiantWorld, StarWorld, LayeredWorld, GasGiantLayeredWorld, BurnManWorld)
GasStarWorldType = Union[GasGiantWorld, StarWorld]
BurnmanWorldType = BurnManWorld
LayeredWorldType = Union[BurnManWorld, LayeredWorld, GasGiantLayeredWorld]
TidalWorldType = Union[GasStarWorldType, BurnmanWorldType, TidalWorld, LayeredWorld]
all_tidal_world_types = (GasGiantWorld, StarWorld, TidalWorld, LayeredWorld, BurnmanWorldType, GasGiantLayeredWorld)

world_types = {
    'star'             : StarWorld,
    'gas'              : GasGiantWorld,
    'gas_giant'        : GasGiantWorld,
    'gas_layered'      : GasGiantLayeredWorld,
    'layered_gas_giant': GasGiantLayeredWorld,
    'layered_ice_giant': GasGiantLayeredWorld,
    'ice_giant'        : GasGiantWorld,
    'simple_tide'      : TidalWorld,
    'simple_tidal'     : TidalWorld,
    'layered'          : LayeredWorld,
    'burnman'          : BurnManWorld
}
