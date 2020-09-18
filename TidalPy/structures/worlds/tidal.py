from typing import Union

from .basic import BaseWorld
from ... import debug_mode
from ...exceptions import ConfigPropertyChangeError, OuterscopePropertySetError, \
    InitiatedPropertyChangeError, InnerscopePropertySetError
from ...rheology.complexCompliance.compliance_models import fixed_q, fixed_q_array
from ...tides import SimpleTides
from ...utilities.types import NoneType


class TidalWorld(BaseWorld):

    """ TidalWorld - Provides a simple base to build tidally dissipative worlds off of.

    Adds basic functionality for CPL or CTL worlds.

    See Also
    --------
    Parent Class:
        TidalPy.structures.worlds.BaseWorld
    Child Classes:
        TidalPy.structures.world.GasGiantWorld
        TidalPy.structures.world.StarWorld
        TidalPy.structures.world.LayeredWorld
        TidalPy.structures.world.GasGiantLayeredWorld
        TidalPy.structures.world.BurnManWorld
    """

    world_class = 'simple_tidal'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):
        """ TidalWorld constructor

        Parameters
        ----------
        world_config : dict
            Configuration file used to build the world. User provided configs override default configurations that
                TidalPy assumes.
            Please see files stored in <TidalPy directory>/structures/worldConfigs for example configuration dict.
        name : str = None
            Name of the world. If None, will use name provided in world_config.
        initialize : bool = True
            Determines if initial reinit should be performed on the world (loading in data from world_config).
        """

        super().__init__(world_config, name=name, initialize=False)

        # Configuration properties
        self.tidal_scale = 1.
        self._is_spin_sync = False
        self._tides_on = False

        # State Properties
        self._fixed_q = None
        self._fixed_dt = None

        # Initialized properties
        self._tides = None  # type: Union[NoneType, SimpleTides]

        # Helper Functions
        # TODO: why are these stored here?
        self.fixed_q_func = fixed_q
        self.fixed_q_func_array = fixed_q_array

        if initialize:
            self.reinit(initial_init=True, setup_simple_tides=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, set_by_burnman: bool = False,
               setup_simple_tides: bool = True):
        """ Initialize or Reinitialize the world based on changes to its configurations.

        This must be called at least once before an instance can be used. The constructor will automatically make an
            initial call to reinit unless told to not to.

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this function has been called.
        reinit_geometry : bool = True
            If `True`, the initializer will automatically call the `set_geometry()` method.
        set_by_burnman : bool = False
            Set to `True` if called from a burnman world.
        setup_simple_tides : bool = True
            Set to `True` if a global CPL/CTL tidal calculation is desired.
        """

        super().reinit(initial_init=initial_init, reinit_geometry=reinit_geometry, set_by_burnman=set_by_burnman)

        # Load in configurations
        self._is_spin_sync = self.config['force_spin_sync']
        self._fixed_dt = self.config['fixed_time_lag']
        self._fixed_q = self.config['quality_factor']
        self._tides_on = self.config['tides_on']

        # Setup the simple tides class - this first switch determines if the world should even check if tides should
        #    be on or not (to allow for child classes to reinit different kinds of tides models)
        if setup_simple_tides:
            if self.tides_on:
                self._tides = SimpleTides(self, store_config_in_world=self.config['store_tides_config_in_world'])
            else:
                self._tides = None

    def clear_state(self, preserve_orbit: bool = False):

        super().clear_state(preserve_orbit)

        self.tides.clear_state()

    def set_fixed_q(self, new_fixed_q: float, run_updates: bool = True):
        """ Set a new global tidal dissipation efficiency or 'fixed-Q' for the world

        Notes
        -----
        .. There is some discussion that, mathematically, q should never be below 1. TidalPy does not check for this
            but it is something to keep in mind.

        Parameters
        ----------
        new_fixed_q : float
            New dissipation efficiency.
        run_updates : bool = True
            If `True`, this method will call the update tides method.
        """

        if debug_mode:
            assert type(new_fixed_q) == float

        self._fixed_q = new_fixed_q

        if run_updates:
            self.update_tides()

    def set_fixed_dt(self, new_fixed_dt: float, run_updates: bool = True):
        """ Set a new global tidal dissipation efficiency frequency scale 'dt' for the world
            (this is used for the CTL model)

        Parameters
        ----------
        new_fixed_dt : float
            New dissipation efficiency frequency scale [s].
        run_updates : bool = True
            If `True`, this method will call the update tides method.
        """

        if new_fixed_dt:
            assert type(new_fixed_dt) == float

        self._fixed_dt = new_fixed_dt

        if run_updates:
            self.update_tides()

    def update_orbit_spin(self, frequency_change: bool = True):
        """ Performs state updates whenever the planet's orbital parameters are changed

        Parameters
        ----------
        frequency_change : bool = True
            If the orbit or spin rate changed then additional calculations must take place.
        """

        super().update_orbit()

        # A change to the orbit will also change tidal frequencies and eccentricity funcs
        self.tides.update_orbit_spin(frequency_change=frequency_change)

        self.update_tides()

    def update_tides(self):
        """ Tell the tides model to calculate tidal heating and potential derivatives """

        super().update_tides()

        self.tides.collapse_modes()


    # # Initialized properties
    @property
    def tides(self):
        """ World's `Tides` class instance - used to calculate all tidal parameters. """
        return self._tides

    @tides.setter
    def tides(self, value):
        raise InitiatedPropertyChangeError


    # # Configuration properties
    @property
    def is_spin_sync(self) -> bool:
        """ If `True`, then the world will be forced in to synchronous rotation """
        return self._is_spin_sync

    @is_spin_sync.setter
    def is_spin_sync(self, value: bool):
        raise ConfigPropertyChangeError

    @property
    def tides_on(self) -> bool:
        """ If `False`, then tides will not be calculated for this world """
        return self._tides_on

    @tides_on.setter
    def tides_on(self, value):
        raise ConfigPropertyChangeError


    # # State properties
    @property
    def fixed_q(self) -> float:
        """ World's global effective tidal dissipation efficiency 'fixed-Q' """
        return self._fixed_q

    @fixed_q.setter
    def fixed_q(self, new_fixed_q):
        self.set_fixed_q(new_fixed_q)

    @property
    def fixed_dt(self) -> float:
        """ World's global tidal dissipation efficiency frequency scale """
        return self._fixed_dt

    @fixed_dt.setter
    def fixed_dt(self, new_fixed_dt):
        self.set_fixed_dt(new_fixed_dt)


    # # Inner-scope properties - Tides model
    @property
    def eccentricity_truncation_lvl(self):
        """ Inner-scope wrapper for tides.eccentricity_truncation_lvl """
        return self.tides.eccentricity_truncation_lvl

    @eccentricity_truncation_lvl.setter
    def eccentricity_truncation_lvl(self, value):
        self.tides.eccentricity_truncation_lvl = value

    @property
    def max_tidal_order_lvl(self):
        """ Inner-scope wrapper for tides.max_tidal_order_lvl """
        return self.tides.max_tidal_order_lvl

    @max_tidal_order_lvl.setter
    def max_tidal_order_lvl(self, value):
        self.tides.max_tidal_order_lvl = value

    @property
    def tidal_susceptibility_reduced(self):
        """ Inner-scope wrapper for tides.tidal_susceptibility_reduced """
        return self.tides.tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        raise InnerscopePropertySetError

    @property
    def tidal_susceptibility(self):
        """ Inner-scope wrapper for tides.tidal_susceptibility """
        return self.tides.tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise InnerscopePropertySetError

    @property
    def unique_tidal_frequencies(self):
        """ Inner-scope wrapper for tides.unique_tidal_frequencies """
        return self.tides.unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise InnerscopePropertySetError

    @property
    def tidal_terms_by_frequency(self):
        """ Inner-scope wrapper for tides.tidal_terms_by_frequency """
        return self.tides.tidal_terms_by_frequency

    @tidal_terms_by_frequency.setter
    def tidal_terms_by_frequency(self, value):
        raise InnerscopePropertySetError

    @property
    def tidal_heating_global(self):
        """ Inner-scope wrapper for tides.tidal_heating_global """
        return self.tides.tidal_heating_global

    @tidal_heating_global.setter
    def tidal_heating_global(self, value):
        raise InnerscopePropertySetError

    @property
    def negative_imk_global(self):
        """ Inner-scope wrapper for tides.negative_imk_global """
        return self.tides.negative_imk_global

    @negative_imk_global.setter
    def negative_imk_global(self, value):
        raise InnerscopePropertySetError

    @property
    def dUdM(self):
        """ Inner-scope wrapper for tides.dUdM """
        return self.tides.dUdM

    @dUdM.setter
    def dUdM(self, value):
        raise InnerscopePropertySetError

    @property
    def dUdw(self):
        """ Inner-scope wrapper for tides.dUdw """
        return self.tides.dUdw

    @dUdw.setter
    def dUdw(self, value):
        raise InnerscopePropertySetError

    @property
    def dUdO(self):
        """ Inner-scope wrapper for tides.dUdO """
        return self.tides.dUdO

    @dUdO.setter
    def dUdO(self, value):
        raise InnerscopePropertySetError

    @property
    def tidal_polar_torque(self):
        """ Inner-scope wrapper for tides.tidal_polar_torque """
        return self.tides.tidal_polar_torque

    @tidal_polar_torque.setter
    def tidal_polar_torque(self, value):
        raise InnerscopePropertySetError

    @property
    def spin_rate_derivative(self):
        """ Inner-scope wrapper for tides.spin_rate_derivative """
        return self.tides.spin_rate_derivative

    @spin_rate_derivative.setter
    def spin_rate_derivative(self, value):
        raise InnerscopePropertySetError


    # # Aliased properties
    @property
    def fixed_time_lag(self):
        """ Alias for self.fixed_dt """
        return self.fixed_dt

    @fixed_time_lag.setter
    def fixed_time_lag(self, new_fixed_time_lag):
        self.fixed_dt = new_fixed_time_lag

    @property
    def use_nsr(self):
        """ Alias for self.is_spin_sync """
        return not self.is_spin_sync

    @use_nsr.setter
    def use_nsr(self, value):
        self.is_spin_sync = ~value

    @property
    def tidal_heating(self):
        """ Alias for self.tidal_heating_global """
        return self.tidal_heating_global

    @tidal_heating.setter
    def tidal_heating(self, value):
        self.tidal_heating_global = value

