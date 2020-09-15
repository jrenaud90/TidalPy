from .layered import LayeredWorld
from .tidal import TidalWorld
from ...exceptions import NotYetImplementedError


# TODO: Implement a fixed-q tides class/method for stellar and gas planets. Wait it is a tidal world...

class GasGiantWorld(TidalWorld):

    world_class = 'gas_giant'


class GasGiantLayeredWorld(LayeredWorld):

    world_class = 'gas_giant_layered'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):

        raise NotYetImplementedError('Layered Gas Giant worlds are not yet implemented. You could try to hack a regular LayeredWorld instead.')
        super().__init__(world_config, name=name, initialize=initialize)