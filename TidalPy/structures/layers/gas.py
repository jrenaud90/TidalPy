from typing import TYPE_CHECKING

from .basic import LayerBase

if TYPE_CHECKING:
    from ..worlds import GasStarWorldType

class GasLayer(LayerBase):

    """ GasLayer class used to store properties and methods for gaseous (and plasma) layers.

    Currently this is unused as tidal heating within gas planets is handled at the world.tides level (since most
        studies assume a global CPL/CTL model for these worlds). This class is a placeholder for future development.
    """

    def __init__(self, layer_name: str, layer_index: int, world: 'GasStarWorldType', layer_config: dict,
                 initialize: bool = True):
        super().__init__(layer_name, layer_index, world, layer_config, initialize)
