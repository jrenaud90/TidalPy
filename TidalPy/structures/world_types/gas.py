from TidalPy.exceptions import NotYetImplementedError

from .layered import LayeredWorld
from .tidal import TidalWorld


class GasGiantWorld(TidalWorld):
    """ GasGiantWorld
    Worlds that are simple gas giants and dissipate tidal energy through the CPL/CTL method (or not at all).


    See Also
    --------
    Parent Class:
        TidalPy.structures.world_types.TidalWorld
    """

    world_class = 'gas_giant'


class GasGiantLayeredWorld(LayeredWorld):
    """ GasGiantLayeredWorld
    Worlds that are gas or ice giants and dissipate tidal energy through either the CPL/CTL method or a more complex
        rheology.

    These world types are not implemented as of at lease v0.2.1


    See Also
    --------
    Parent Class:
        TidalPy.structures.world_types.LayeredWorld
    """

    world_class = 'gas_giant_layered'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):
        raise NotYetImplementedError(
            'Layered Gas Giant world_types are not yet implemented. You could try to hack a regular LayeredWorld instead.'
            )
