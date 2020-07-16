from exceptions import IncorrectAttributeType, ConfigPropertyChangeError, OuterscopePropertySetError
from .basic import BaseWorld
from ...rheology.complexCompliance.compliance_models import fixed_q, fixed_q_array
from ...tides.tides import SimpleTides


class TidalWorld(BaseWorld):

    """ TidalWorld - Provides a simple base to build tidally dissipative worlds off of.

    Adds basic functionality for CPL or CTL worlds.
    """

    world_class = 'simple_tidal_world'

    def __init__(self, planet_config: dict, name: str = None, initialize: bool = True):

        super().__init__(planet_config, name=name, initialize=False)

        # State Properties
        self._fixed_q = None
        self._fixed_dt = None
        self._tidal_host = None

        # Configurations
        self._is_spin_sync = False

        # Class Objects
        self.tides = None  # type: SimpleTides

        # Helper Functions
        self.fixed_q_func = fixed_q
        self.fixed_q_func_array = fixed_q_array

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):

        super().reinit(initial_init)

        # Load in configurations
        self._is_spin_sync = self.config['is_spin_sync']
        self._fixed_dt = self.config['fixed_time_lag']
        self._fixed_q = self.config['quality_factor']

        # Setup the simple tides class
        self.tides = SimpleTides(self, store_config_in_world=self.config['store_tides_config_in_world'])


    # State properties
    @property
    def fixed_q(self) -> float:
        return self._fixed_q

    @fixed_q.setter
    def fixed_q(self, new_fixed_q: float):

        if type(new_fixed_q) is not float:
            raise IncorrectAttributeType

        self._fixed_q = new_fixed_q
        self.update_tides()

    @property
    def fixed_dt(self) -> float:
        return self._fixed_dt

    @fixed_dt.setter
    def fixed_dt(self, new_fixed_time_lag: float):

        if type(new_fixed_time_lag) is not float:
            raise IncorrectAttributeType

        self._fixed_dt = new_fixed_time_lag
        self.update_tides()

    @property
    def is_spin_sync(self) -> bool:
        return self._is_spin_sync

    @is_spin_sync.setter
    def is_spin_sync(self, value: bool):
        raise ConfigPropertyChangeError

    # Inner-scope properties - Tides model
    @property
    def eccentricity_truncation_lvl(self):
        return self.tides.eccentricity_truncation_lvl

    @eccentricity_truncation_lvl.setter
    def eccentricity_truncation_lvl(self, value):
        self.tides.eccentricity_truncation_lvl = value

    @property
    def max_tidal_order_lvl(self):
        return self.tides.max_tidal_order_lvl

    @max_tidal_order_lvl.setter
    def max_tidal_order_lvl(self, value):
        self.tides.max_tidal_order_lvl = value

    @property
    def tidal_susceptibility_reduced(self):
        return self.tides.tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        self.tides.tidal_susceptibility_reduced = value

    @property
    def tidal_susceptibility(self):
        return self.tides.tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        self.tides.tidal_susceptibility = value

    @property
    def unique_tidal_frequencies(self):
        return self.tides.unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        self.tides.unique_tidal_frequencies = value

    @property
    def tidal_terms_by_frequency(self):
        return self.tides.tidal_terms_by_frequency

    @tidal_terms_by_frequency.setter
    def tidal_terms_by_frequency(self, value):
        self.tides.tidal_terms_by_frequency = value

    @property
    def tidal_heating_global(self):
        return self.tides.tidal_heating_global

    @tidal_heating_global.setter
    def tidal_heating_global(self, value):
        self.tides.tidal_heating_global = value

    @property
    def negative_imk_global(self):
        return self.tides.negative_imk_global

    @negative_imk_global.setter
    def negative_imk_global(self, value):
        self.tides.negative_imk_global = value

    @property
    def dUdM(self):
        return self.tides.dUdM

    @dUdM.setter
    def dUdM(self, value):
        self.tides.dUdM = value

    @property
    def dUdw(self):
        return self.tides.dUdw

    @dUdw.setter
    def dUdw(self, value):
        self.tides.dUdw = value

    @property
    def dUdO(self):
        return self.tides.dUdO

    @dUdO.setter
    def dUdO(self, value):
        self.tides.dUdO = value

    @property
    def tidal_polar_torque(self):
        return self.tides.tidal_polar_torque

    @tidal_polar_torque.setter
    def tidal_polar_torque(self, value):
        self.tides.tidal_polar_torque = value

    @property
    def spin_rate_derivative(self):
        return self.tides.spin_rate_derivative

    @spin_rate_derivative.setter
    def spin_rate_derivative(self, value):
        self.tides.spin_rate_derivative = value

    # Outer-scope properties
    @property
    def tidal_host(self):
        if self.orbit is None:
            return None
        else:
            return self.tidal_hosts[self]

    @tidal_host.setter
    def tidal_host(self, value):
        raise OuterscopePropertySetError


    # Alias properties
    @property
    def fixed_time_lag(self):
        return self.fixed_dt

    @fixed_time_lag.setter
    def fixed_time_lag(self, new_fixed_time_lag):
        self.fixed_dt = new_fixed_time_lag

    @property
    def use_nsr(self):
        return not self.is_spin_sync

    @use_nsr.setter
    def use_nsr(self, value):
        self.is_spin_sync = ~value

    @property
    def tidal_heating(self):
        return self.tidal_heating_global

    @tidal_heating.setter
    def tidal_heating(self, value):
        self.tidal_heating_global = value

