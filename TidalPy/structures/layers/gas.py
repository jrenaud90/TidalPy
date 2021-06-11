from typing import TYPE_CHECKING

from .basic import LayerBase

if TYPE_CHECKING:
    from ..world_types import GasGiantLayeredWorld


class GasLayer(LayerBase):
    """" GasLayer
    Layer object used to construct gas giant or ice giant planets that contain a significant gas layer. Currently,
        these layers do not do anything over the base layer class but allow for future functionality.

    Notes:
    .. Does not provide any functionality to perform tidal calculations (see PhysicsLayer instead)

    See Also
    --------
    TidalPy.structures.layers.LayerBase
    """

    layer_class = 'gas'

    def __init__(
        self, layer_name: str, layer_index: int, world: 'GasGiantLayeredWorld', layer_config: dict,
        is_top_layer: bool, initialize: bool = True
        ):
        """ Gas layer constructor

        Parameters
        ----------
        layer_name : str
            User-friendly name of layer.
        layer_index : int
            Location of layer within a world (0 indicates center-most).
        world : LayeredWorldType
            World instance where layer was initialized in.
        layer_config : dict
            Layer's user-provided configurations.
        is_top_layer : bool
            If `True`, this layer is the top-most layer.
        initialize : bool = True
            If `True`, then the Layer's reinit is called at the end of the constructor.
        """

        super().__init__(layer_name, layer_index, world, layer_config, is_top_layer, initialize)
