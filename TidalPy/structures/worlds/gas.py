from .tidal import TidalWorld

# TODO: Implement a fixed-q tides class/method for stellar and gas planets. Wait it is a tidal world...

class GasGiantWorld(TidalWorld):

    world_class = 'gas_giant'
