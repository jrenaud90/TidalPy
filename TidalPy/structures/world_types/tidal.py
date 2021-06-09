from typing import Union

from .basic import BaseWorld
from ... import log
from ...exceptions import ConfigPropertyChangeError, InitiatedPropertyChangeError, InnerscopePropertySetError, \
    AttributeNotSetError, IncorrectMethodToSetStateProperty
from ...rheology.complex_compliance.compliance_models import fixed_q
from ...tides import GlobalApproxTides
from ...utilities.types import NoneType, FloatArray


class TidalWorld(BaseWorld):

    """ TidalWorld - Provides a simple base to build tidally dissipative world_types off of.

    Adds basic functionality for CPL or CTL world_types.

    See Also
    --------
    Parent Class:
        TidalPy.structures.world_types.BaseWorld
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
            Please see files stored in <TidalPy directory>/structures/world_configs for example configuration dict.
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
        self._spin_time_derivative = None
        self._tidal_polar_torque = None

        # Initialized properties
        self._tides = None  # type: Union[NoneType, GlobalApproxTides]

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
        self._tides_on = self.config['tides_on']

        # Setup the simple tides class - this first switch determines if the world should even check if tides should
        #    be on or not (to allow for child methods to reinit different kinds of tides models)
        if setup_simple_tides:
            if self.tides_on:
                self._tides = GlobalApproxTides(self, store_config_in_world=self.config['store_tides_config_in_world'])
            else:
                self._tides = None

    def orbit_spin_changed(self, orbital_freq_changed: bool = False, spin_freq_changed: bool = False,
                           eccentricity_changed: bool = False, obliquity_changed: bool = False,
                           call_orbit_dissipation: bool = True):
        """ The world's orbit, spin, and/or obliquity has changed. Make any necessary updates """

        super().orbit_spin_changed(orbital_freq_changed=orbital_freq_changed, spin_freq_changed=spin_freq_changed,
                                   eccentricity_changed=eccentricity_changed, obliquity_changed=obliquity_changed,
                                   call_orbit_dissipation=call_orbit_dissipation)

        # Tell the tides model that the orbit and/or spin has changed.
        if self.tides is not None:
            self.tides.orbit_spin_changed(eccentricity_change=eccentricity_changed,
                                          obliquity_change=obliquity_changed,
                                          orbital_freq_changed=orbital_freq_changed,
                                          spin_freq_changed=spin_freq_changed)

            if call_orbit_dissipation:
                self.orbit.dissipation_changed(self)

    def tidal_frequencies_changed(self, collapse_tidal_modes: bool = True):
        """ The tidal frequencies have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        log.debug(f'Tidal frequencies changed for {self}.')

        # Updates set by child methods.

    def dissipation_changed(self):
        """ Tidal dissipation has changed. Make any necessary updates. """

        log.debug(f'Tidal dissipation changed for {self}.')

        if self.orbit is not None:
            self.orbit.dissipation_changed(self)

        self.internal_heating_changed()

    def internal_heating_changed(self):
        """ The internal heating of this world has changed. Make any necessary updates. """

        log.debug(f'Internal heating changed for {self}.')

        if self.world_class not in ['layered', 'burnman']:
            # Layered world_types use a different method to change the surface temperature.
            self.update_surface_temperature()

    def clear_state(self, preserve_orbit: bool = False):
        """ Clear the world's current state variables back to their defaults.

        The defaults may be Nones or set by the user-provided configuration.

        Parameters
        ----------
        preserve_orbit: bool = False
            If `True`, data about this planet's orbit will be cleared from any associated Orbit methods.
        """

        super().clear_state(preserve_orbit)

        self.tides.clear_state()

    def set_fixed_q(self, fixed_q: float, run_updates: bool = True):
        """ Set a new global tidal dissipation efficiency or 'fixed-Q' for the world

        This is used in calculating tidal dissipation assuming a CPL model.

        Notes
        -----
        .. There is some discussion that, mathematically, q should never be below 1. TidalPy does not check for this
            but it is something to keep in mind.

        Parameters
        ----------
        fixed_q : float
            New dissipation efficiency.
        run_updates : bool = True
            If `True`, this method will call the update tides method.
        """

        if self.tides is not None:
            self.tides.set_fixed_q(fixed_q, run_updates=run_updates)

    def set_fixed_dt(self, fixed_dt: float, run_updates: bool = True):
        """ Set a new global tidal dissipation efficiency frequency scale 'dt' for the world

        This is used in calculating tidal dissipation assuming a CTL model.

        Parameters
        ----------
        fixed_dt : float
            New dissipation efficiency frequency scale [s].
        run_updates : bool = True
            If `True`, this method will call the update tides method.
        """

        if self.tides is not None:
            self.tides.set_fixed_dt(fixed_dt, run_updates=run_updates)

    def calc_spin_derivative(self) -> FloatArray:
        """ Calculate spin-rate time derivative based on the world's current state.

        Requires the tides model to have already calculated dU / dO
            (tidal potential derivative with respect to the node)

        Returns
        -------
        spin_rate_derivative : FloatArray
            Spin-rate derivative for the planet [rads s-2]
        """

        if self.orbit is not None:

            if self.dUdO is None:
                raise AttributeNotSetError(f'Potential derivatives not calculated for {self}.')

            self._spin_time_derivative = self.tidal_host.mass * self.dUdO / self.moi
            self._tidal_polar_torque = self.tidal_host.mass * self.dUdO

        return self.spin_time_derivative

    def get_internal_heating_to_surface(self) -> Union[NoneType, FloatArray]:
        """ Get the amount of internal heating that is making it to the surface.

        Returns
        -------
        internal_heating_to_surface : Union[NoneType, FloatArray]
            Amount of heating that is reaching the surface [W]
        """

        internal_heating_to_surface = None
        if self.tides is not None:
            if self.tides.tidal_heating_global is not None:
                internal_heating_to_surface = self.tides.tidal_heating_global * self.internal_to_surf_heating_frac

        return internal_heating_to_surface


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
    def spin_time_derivative(self) -> FloatArray:
        """ The time derivative of the spin rate [rad s-2] """
        return self._spin_time_derivative

    @spin_time_derivative.setter
    def spin_time_derivative(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def tidal_polar_torque(self) -> FloatArray:
        """ The tidal polar torque [N m] """
        return self._tidal_polar_torque

    @tidal_polar_torque.setter
    def tidal_polar_torque(self, value):
        raise IncorrectMethodToSetStateProperty


    # # Inner-scope properties - Tides model
    @property
    def fixed_q(self) -> Union[NoneType, float]:
        """ World's global effective tidal dissipation efficiency 'fixed-Q' """
        if self.tides is None:
            return None
        else:
            return self.tides.fixed_q

    @fixed_q.setter
    def fixed_q(self, new_fixed_q):
        self.tides.set_fixed_q(new_fixed_q)

    @property
    def fixed_dt(self) -> Union[NoneType, float]:
        """ World's global tidal dissipation efficiency frequency scale """
        if self.tides is None:
            return None
        else:
            return self.tides.fixed_dt

    @fixed_dt.setter
    def fixed_dt(self, new_fixed_dt):
        self.tides.set_fixed_dt(new_fixed_dt)

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
    def global_negative_imk_by_orderl(self):
        """ Global negative of the imaginary portion of the Love number, -Im[k_l] """
        return self.tides.global_negative_imk_by_orderl

    @global_negative_imk_by_orderl.setter
    def global_negative_imk_by_orderl(self, value):
        raise InnerscopePropertySetError

    @property
    def global_love_by_orderl(self):
        """ Global complex Love number, k_l """
        return self.tides.global_love_by_orderl

    @global_love_by_orderl.setter
    def global_love_by_orderl(self, value):
        raise InnerscopePropertySetError

    @property
    def effective_q_by_orderl(self):
        """ World's effective tidal dissipation factor (for each tidal order level) """
        return self.tides.effective_q_by_orderl

    @effective_q_by_orderl.setter
    def effective_q_by_orderl(self, value):
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

