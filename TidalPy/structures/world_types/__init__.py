from typing import Union

from .basic import BaseWorld
from .gas import GasGiantLayeredWorld, GasGiantWorld
from .layered import LayeredWorld
from .stellar import StarWorld
from .tidal import TidalWorld

AllWorldType = Union[BaseWorld, TidalWorld, GasGiantWorld, StarWorld, LayeredWorld, GasGiantLayeredWorld]
GasStarWorldType = Union[GasGiantWorld, StarWorld]
LayeredWorldType = Union[LayeredWorld, GasGiantLayeredWorld]
TidalWorldType = Union[GasStarWorldType, TidalWorld, LayeredWorld]

all_world_types = (BaseWorld, TidalWorld, GasGiantWorld, StarWorld, LayeredWorld, GasGiantLayeredWorld)
all_tidal_world_types = (GasGiantWorld, StarWorld, TidalWorld, LayeredWorld, GasGiantLayeredWorld)

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
    'layered'          : LayeredWorld
    }
