from .tidal import TidalWorld

class ThermalWorld(TidalWorld):

    world_class = 'thermal'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):

        super().__init__(world_config, name=name, initialize=False)

        if initialize:
            self.reinit()