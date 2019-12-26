from . import GeometricWorld

from ...tides.tides import Tides


class TidalWorld(GeometricWorld):

    world_class = 'tidal'

    def __init__(self, world_config: dict, name: str = None):

        super().__init__(world_config, name=name)

        # Things setup in reinit
        self.tides = None

    def reinit(self):

        super().reinit()

        # Setup tidal model
        self.tides = Tides(self)


    # Inner scope properties - Tides model
    @property
    def tidal_modes(self):
        return self.tides.modes

    @tidal_modes.setter
    def tidal_modes(self, value):
        self.tides.modes = value

    @property
    def tidal_frequencies(self):
        return self.tides.frequencies

    @tidal_frequencies.setter
    def tidal_frequencies(self, value):
        self.tides.frequencies = value

    @property
    def tidal_unique_tidal_frequencies(self):
        return self.tides.unique_tidal_frequencies

    @tidal_unique_tidal_frequencies.setter
    def tidal_unique_tidal_frequencies(self, value):
        self.tides.unique_tidal_frequencies = value

    @property
    def tidal_mode_names(self):
        return self.tides.mode_names

    @tidal_mode_names.setter
    def tidal_mode_names(self, value):
        self.tides.mode_names = value

    @property
    def tidal_heating_subterms(self):
        return self.tides.heating_subterms

    @tidal_heating_subterms.setter
    def tidal_heating_subterms(self, value):
        self.tides.heating_subterms = value

    @property
    def tidal_ztorque_subterms(self):
        return self.tides.ztorque_subterms

    @tidal_ztorque_subterms.setter
    def tidal_ztorque_subterms(self, value):
        self.tides.ztorque_subterms = value