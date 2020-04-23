from . import GeometricWorld
from ...exceptions import ImproperPropertyHandling, ConfigPropertyChangeError

from ...tides.tides import Tides


class SimpleTidalWorld(GeometricWorld):
    """ TidalWorld - Provides a simple base to build tidally dissipative worlds off of.
    """

    world_class = 'simple_tide'

    def __init__(self, planet_config: dict, name: str = None, initialize: bool = True):
        super().__init__(planet_config, name=name, initialize=False)

        # State Properties

        # Additional models to be set in reinit
        self.tides = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        super().reinit(initial_init)

        self.tides = Tides(self, store_config_in_world=True)


class TidalWorld(SimpleTidalWorld):

    world_class = 'tidal'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):

        super().__init__(world_config, name=name, initialize=False)

        # State properties
        self._is_spin_sync = None

        # Things setup in reinit
        self._tides = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):

        super().reinit(initial_init)

        # Load in configurations
        self._is_spin_sync = self.config['is_spin_sync']

        # Setup tidal model
        self._tides = Tides(self)

    # State properties
    @property
    def tides(self):
        return self._tides

    @tides.setter
    def tides(self, value):
        raise ConfigPropertyChangeError

    @property
    def is_spin_sync(self):
        return self._is_spin_sync

    @is_spin_sync.setter
    def is_spin_sync(self, value):
        raise ConfigPropertyChangeError

    # Inner scope properties - Tides model
    @property
    def orbital_truncation_lvl(self):
        return self.tides.orbital_truncation_lvl

    @orbital_truncation_lvl.setter
    def orbital_truncation_lvl(self, value):
        self.tides.orbital_truncation_lvl = value

    @property
    def tidal_order_lvl(self):
        return self.tides.tidal_order_lvl

    @tidal_order_lvl.setter
    def tidal_order_lvl(self, value):
        self.tides.tidal_order_lvl = value

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

    @property
    def tidal_heating_bylayer(self):
        return self.tides.tidal_heating_bylayer

    @tidal_heating_bylayer.setter
    def tidal_heating_bylayer(self, value):
        self.tides.tidal_heating_bylayer = value

    @property
    def tidal_ztorque_bylayer(self):
        return self.tides.tidal_ztorque_bylayer

    @tidal_ztorque_bylayer.setter
    def tidal_ztorque_bylayer(self, value):
        self.tides.tidal_ztorque_bylayer = value

    # Alias properties
    @property
    def use_nsr(self):
        return not self.is_spin_sync

    @use_nsr.setter
    def use_nsr(self, value):
        raise ImproperPropertyHandling

