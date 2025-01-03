""" Base Tides Module
"""

import numpy as np
from typing import Dict, TYPE_CHECKING, Tuple

import TidalPy
from TidalPy.exceptions import (AttributeNotSetError, ConfigPropertyChangeError, IncompatibleModelError,
                                IncorrectMethodToSetStateProperty, InitiatedPropertyChangeError, NotYetImplementedError,
                                OuterscopePropertySetError, ParameterValueError)
from TidalPy.utilities.classes.config.config import WorldConfigHolder


from ..love1d import complex_love_general, effective_rigidity_general
from ..dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced
from ..modes.mode_manipulation import find_mode_manipulators

from TidalPy.logger import get_logger
log = get_logger("TidalPy")


if TYPE_CHECKING:
    from TidalPy.utilities.types import ComplexArray, FloatArray
    from TidalPy.structures.world_types import TidalWorldType

    from ..eccentricity_funcs import EccenOutput
    from ..inclination_funcs import InclinOutput
    from ..modes.mode_manipulation import (DissipTermsArray, FreqSig, ResultsByFreqType, UniqueFreqType)


# TODO: Add a spin-sync version


class TidesBase(WorldConfigHolder):
    """ TidesBase
    Holder for all tidal heating and tidal potential calculations

    Tides class stores model parameters and methods for calculating tidal heating and tidal potential derivatives
        which are general functions of (T, P, melt_frac, w, e, obliquity)

    Attributes
    ----------
    eccentricity_truncation_lvl
    max_tidal_order_lvl
    use_obliquity_tides
    multiply_modes_by_sign
    eccentricity_results
    obliquity_results
    tidal_susceptibility_reduced
    tidal_susceptibility
    unique_tidal_frequencies
    tidal_terms_by_frequency
    tidal_heating_global
    global_negative_imk_by_orderl
    global_love_by_orderl
    effective_q_by_orderl
    dUdM
    dUdw
    dUdO
    fixed_q
    fixed_dt
    fixed_k2
    thermal_set
    orbit_set

    """

    model = 'base'
    world_config_key = 'tides'

    def __init__(self, world: 'TidalWorldType', store_config_in_world: bool = True, initialize: bool = True):
        """ Constructor for TidesBase class

        Parameters
        ----------
        world : TidalWorldType
            The world where tides are being calculated.
        store_config_in_world : bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `world.config` dictionary.
        initialize : bool = True
            If `True`, then an initial call to the tide's reinit method will be made at the end of construction.
        """

        self.default_config = TidalPy.config['tides']['models'][self.model]

        super().__init__(world, store_config_in_world=store_config_in_world)

        # Initialized properties
        self._thermal_set = False
        self._orbit_set = False

        # State properties
        self._fixed_q = None
        self._fixed_dt = None
        self._fixed_k2 = None
        self._tidal_susceptibility = None
        self._tidal_susceptibility_reduced = None
        self._unique_tidal_frequencies = None
        self._tidal_terms_by_frequency = None
        self._tidal_heating_global = None
        self._dUdM = None
        self._dUdw = None
        self._dUdO = None
        self._effective_q_by_orderl = None
        self._global_negative_imk_by_orderl = None
        self._global_love_by_orderl = None
        self._need_to_collapse_modes = False
        self._new_tidal_frequencies = False

        # Model configurations that will be set in reinit
        self._eccentricity_truncation_lvl = None
        self._max_tidal_order_lvl = None
        self._use_obliquity_tides = None
        self._multiply_modes_by_sign = None

        # Functions to be initialized in self.reinit
        self._eccentricity_results = None
        self._obliquity_results = None
        self.eccentricity_func = None
        self.obliquity_func = None
        self.collapse_modes_func = None
        self.calculate_modes_func = None

        # TidalPy logging and debug info
        log.debug(f'Building {self.model} tides class for {self.world.name}.')

        # Call reinit for initialization
        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Load configurations into the Tides class and import any config-dependent functions.

        This reinit process is separate from the __init__ method because the Orbit class may need to overload some
            configurations after class initialization.

        Parameters
        ----------
        initial_init : bool = False
            This should be set to True the first time reinit is called.

        Raises
        ------
        NotYetImplementedError
        ParameterValueError

        """

        if initial_init:
            log.debug(f'Initializing tides class: {self}.')
        else:
            log.debug(f'Reinit called for {self}.')
            self.clear_state()

            # Reset configurations
            self._eccentricity_truncation_lvl = None
            self._max_tidal_order_lvl = None
            self._use_obliquity_tides = None
            self._multiply_modes_by_sign = None

        # Load in configurations
        self._eccentricity_truncation_lvl = self.config['eccentricity_truncation_lvl']
        self._max_tidal_order_lvl         = self.config['max_tidal_order_l']
        self._use_obliquity_tides         = self.config['obliquity_tides_on']
        self._multiply_modes_by_sign      = self.config['multiply_modes_by_sign']

        # Ensure the tidal order and orbital truncation levels make sense
        if self.max_tidal_order_lvl > 7:
            raise NotYetImplementedError(f'Tidal order {self.max_tidal_order_lvl} has not been implemented.')
        if self.eccentricity_truncation_lvl % 2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self.eccentricity_truncation_lvl < 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self.eccentricity_truncation_lvl not in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
            if self.max_tidal_order_lvl == 2 and self.eccentricity_truncation_lvl == 22:
                # This was implemented in v0.2.1
                pass
            else:
                raise NotYetImplementedError(
                    f'Orbital truncation level of {self.eccentricity_truncation_lvl} has not been implemented.'
                    )

        # Setup functions
        self.calculate_modes_func, self.collapse_modes_func, self.eccentricity_func, self.obliquity_func = \
            find_mode_manipulators(self.max_tidal_order_lvl, self.eccentricity_truncation_lvl, self.use_obliquity_tides)

    def post_orbit_initialize(self):
        """ Initialize various tidal parameters once a tidal host is connected to the world (via an orbit class).

        Raises
        ------
        AttributeNotSetError

        """

        log.debug(f'Post-orbit tidal initialization started for {self}.')

        # Pull out parameters that are called multiple times in this method
        tidal_host = self.tidal_host
        world_radius = self.world.radius
        semi_major_axis = self.semi_major_axis

        if tidal_host is None:
            raise AttributeNotSetError('Tidal host must be connected to target body in order to initialize tides.')

        # Try to update properties that depend on the orbit and tidal host
        self._tidal_susceptibility_reduced = \
            calc_tidal_susceptibility_reduced(tidal_host.mass, world_radius)

        if semi_major_axis is not None:
            self._tidal_susceptibility = calc_tidal_susceptibility(
                tidal_host.mass, world_radius, semi_major_axis
                )

    def orbit_spin_changed(
        self, eccentricity_change: bool = True, obliquity_change: bool = True,
        orbital_freq_changed: bool = True, spin_freq_changed: bool = True,
        force_obliquity_update: bool = False,
        call_world_frequency_changed: bool = True,
        call_collapse_modes: bool = True
        ) -> Tuple[ 'UniqueFreqType', 'ResultsByFreqType']:
        """ Calculate tidal heating and potential derivative terms based on the current orbital state.

        This will also calculate new unique tidal frequencies which must then be digested by the rheological model
            at each planetary layer.


        Parameters
        ----------
        eccentricity_change : bool = True
            If there was no change in eccentricity (or if the orbit set the eccentricity) set this to False for a
            performance boost. If False, eccentricity functions won't be called.
        obliquity_change : bool = True
            If there was no change in obliquity set this to False for a performance boost.
            If False, obliquity functions won't be called.
        orbital_freq_changed : bool = True
            If there was no change in the orbital frequency set this to False for a performance boost.
        spin_freq_changed : bool = True
            If there was no change in the spin frequency set this to False.
        force_obliquity_update : bool = False
            If True, then this method will call the obliquity update function even if it would otherwise skip it
            due to obliquity dependence being turned off.
        call_world_frequency_changed : bool = True
            If True, then the method will call the world's complex compliance calculator.
            This flag is set to False for the CPL method.
        call_collapse_modes : bool = True
            If True, then this method will call collapse modes (if needed).

        Returns
        -------
        unique_tidal_frequencies : UniqueFreqType
            Each unique frequency stored as a signature (orbital motion and spin-rate coeffs),
            and the calculated frequency (combo of orbital motion and spin-rate) [rad s-1]
        tidal_terms_by_frequency : ResultsByFreqType
            Results for tidal heating, dU/dM, dU/dw, dU/dO are stored in a tuple for each tidal harmonic l and
            unique frequency.

        See Also
        --------
        TidalPy.tides.dissipation.mode_collapse

        """

        log.debug(f'Method orbit_spin_changed called for {self}.')

        self._need_to_collapse_modes = False
        self._new_tidal_frequencies = False

        # Pull out parameters that are called multiple times in the method.
        eccentricity = self.eccentricity
        obliquity = self.obliquity
        world_radius = self.world.radius
        semi_major_axis = self.semi_major_axis
        spin_frequency = self.spin_frequency
        orbital_frequency = self.orbital_frequency

        # Make updates
        if orbital_freq_changed:
            # Orbital changes may have changed the tidal susceptibility
            self._tidal_susceptibility = calc_tidal_susceptibility(
                self.tidal_host.mass, world_radius,
                semi_major_axis
                )

        # Check if we need to calculate new eccentricity results
        if eccentricity is not None:
            if eccentricity_change or self.eccentricity_results is None:
                self._eccentricity_results = self.eccentricity_func(eccentricity)
                self._need_to_collapse_modes = True


        # Check if we need to update obliquity results
        if self.use_obliquity_tides:
            if obliquity is not None:
                if obliquity_change or self.obliquity_results is None:
                    self._obliquity_results = self.obliquity_func(obliquity)
                    self._need_to_collapse_modes = True
        else:
            # We assume that the obliquity has the same shape of eccentricity even when it is not being used.
            # OPT: Check to see if we can universally get away with having this just be a scalar rather than a
            #  potentially large array.
            if eccentricity is not None:
                zero_obliquity = np.zeros_like(eccentricity)

                # Since obliquity tides are not used they will never change. So do not call the function unless forced
                #  to do so.
                if self.obliquity_results is None or force_obliquity_update:
                    self._obliquity_results = self.obliquity_func(zero_obliquity)
                    self._need_to_collapse_modes = True

        obliquity_results = self.obliquity_results
        eccentricity_results = self.eccentricity_results

        # Check the shape of the eccentricity and obliquity results
        if obliquity_results is not None and eccentricity_results is not None:
            # Check that obliquity and eccentricity results have the same size (same max_l used).
            if len(obliquity_results) != len(eccentricity_results):
                raise IncompatibleModelError(
                    'Obliquity and Eccentricity results do not have the same size.'
                    'max_tidal_l may not be equal for both functions.'
                    )

        # Determine if new tidal frequencies need to be calculated
        if spin_freq_changed or orbital_freq_changed:
            if eccentricity_results is not None and obliquity_results is not None and \
                    spin_frequency is not None and orbital_frequency is not None:
                # Update the tidal frequencies and terms using the new orbital frequency
                self._unique_tidal_frequencies, self._tidal_terms_by_frequency = \
                    self.calculate_modes_func(
                        spin_frequency, orbital_frequency, semi_major_axis,
                        world_radius, eccentricity_results, obliquity_results,
                        self.multiply_modes_by_sign
                        )
                # Now that there are new frequencies, tell the world so that new complex compliances can be calculated.
                self._new_tidal_frequencies = True
                self._need_to_collapse_modes = True

        # Tell the planet to make any updates due to a potential frequency change.
        if self._new_tidal_frequencies and call_world_frequency_changed:
            # Updating the tidal frequencies will automatically call the tide's collapse_modes method. The
            #    collapse_tidal_modes flag prevents that method from being called more than once.
            self.world.tidal_frequencies_changed(collapse_tidal_modes=False)

        # Collapse tidal modes
        if self._need_to_collapse_modes and call_collapse_modes:
            self.collapse_modes()

        # Return frequencies and tidal terms
        return self.unique_tidal_frequencies, self.tidal_terms_by_frequency

    def fixed_q_dt_changed(self):
        """ The fixed tidal dissipation parameters (fixed-q or fixed-dt) have changed. Make any necessary updates. """

        log.debug(f'Method fixed_q_dt_changed called for {self}.')

        # Updates done by child classes

    def clear_state(self):
        """ Clear the tidal model's state properties """

        super().clear_state()

        self._eccentricity_results          = None
        self._obliquity_results             = None
        self._tidal_susceptibility          = None
        self._tidal_susceptibility_reduced  = None
        self._unique_tidal_frequencies      = None
        self._tidal_terms_by_frequency      = None
        self._tidal_heating_global          = None
        self._global_negative_imk_by_orderl = None
        self._dUdM = None
        self._dUdw = None
        self._dUdO = None
        self._need_to_collapse_modes = False
        self._new_tidal_frequencies = False

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

        log.debug(f'Method set_fixed_q called for {self}.')

        self._fixed_q = fixed_q

        if run_updates:
            self.fixed_q_dt_changed()

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

        log.debug(f'Method set_fixed_dt called for {self}.')

        self._fixed_dt = fixed_dt

        if run_updates:
            self.fixed_q_dt_changed()

    def set_state(self, fixed_q: float = None, fixed_dt: float = None, run_updates: bool = True):
        """ Change the state of the Tides class

        This is more efficient than calling the individual setters.

        Parameters
        ----------
        fixed_q : float (optional)
            New dissipation efficiency.
        fixed_dt : float (optional)
            New dissipation efficiency frequency scale [s].
        run_updates : bool = True
            If `True`, this method will call the update tides method.
        """

        log.debug(f'Method set_state called for {self}.')

        if fixed_q is not None:
            self.set_fixed_q(fixed_q, run_updates=False)

        if fixed_dt is not None:
            self.set_fixed_dt(fixed_dt, run_updates=False)

        if run_updates:
            self.fixed_q_dt_changed()

    def collapse_modes(self):
        """ Calculate Global Love number based on current thermal state.

        Requires a prior orbit_spin_changed() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        Returns
        -------
        tidal_heating : np.ndarray
            Tidal heating [W]
            This could potentially restricted to a layer or for an entire planet.
        dUdM : np.ndarray
            Tidal potential derivative with respect to the mean anomaly [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.
        dUdw : np.ndarray
            Tidal potential derivative with respect to the pericentre [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.
        dUdO : np.ndarray
            Tidal potential derivative with respect to the planet's node [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.

        See Also
        --------
        TidalPy.tides.Tides.orbit_spin_changed
        """

        log.debug(f'Method collapse_modes called for {self}.')

        # Most of mode_collapse is implemented by a child class of TidesBase

    # # Initialized properties
    @property
    def thermal_set(self) -> bool:
        """ Flag to check if a Thermal Class has been initialized in the tide class host world. """
        return self._thermal_set

    @thermal_set.setter
    def thermal_set(self, value):
        raise InitiatedPropertyChangeError

    @property
    def orbit_set(self) -> bool:
        """ Flag to check if an Orbit Class has been initialized in the tide class host world. """
        return self._orbit_set

    @orbit_set.setter
    def orbit_set(self, value):
        raise InitiatedPropertyChangeError

    # # Configuration properties
    @property
    def eccentricity_truncation_lvl(self) -> int:
        """ Maximum eccentricity truncation level to include in tidal calculations """
        return self._eccentricity_truncation_lvl

    @eccentricity_truncation_lvl.setter
    def eccentricity_truncation_lvl(self, value):
        # TODO: Think about if you want the user to update these. These setters could make a call to self.reinit()
        #    which the user could make on their own. So it may make sense to allow the setter.
        raise ConfigPropertyChangeError

    @property
    def max_tidal_order_lvl(self) -> int:
        """ Maximum tidal order to include in tidal calculations """
        return self._max_tidal_order_lvl

    @max_tidal_order_lvl.setter
    def max_tidal_order_lvl(self, value):
        raise ConfigPropertyChangeError

    @property
    def use_obliquity_tides(self) -> bool:
        """ Flag for if obliquity tides should be calculated

        Notes
        -----
        - Setting this to False leads to more efficient calculations than if the obliquity of a world is simply set to
        zero.
        """

        return self._use_obliquity_tides

    @use_obliquity_tides.setter
    def use_obliquity_tides(self, value):
        raise ConfigPropertyChangeError

    @property
    def multiply_modes_by_sign(self):
        """ Flag for if the tidal dissipation results should be multiplied by the mode's sign.

        This should always be `True` unless doing specific tests.
        """
        return self._multiply_modes_by_sign

    @multiply_modes_by_sign.setter
    def multiply_modes_by_sign(self, value):
        raise ConfigPropertyChangeError

    # # State properties
    @property
    def eccentricity_results(self) -> Dict[int, 'EccenOutput']:
        """ Eccentricity function results (squared) stored by order_l """
        return self._eccentricity_results

    @eccentricity_results.setter
    def eccentricity_results(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def obliquity_results(self) -> Dict[Tuple[int, int], 'InclinOutput']:
        """ Obliquity/Inclination function results (squared) stored by integers (m, p) """
        return self._obliquity_results

    @obliquity_results.setter
    def obliquity_results(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def tidal_susceptibility_reduced(self) -> 'FloatArray':
        """ Tidal susceptibility (reduced, no semi-major axis dependence) [N m7] """
        return self._tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def tidal_susceptibility(self) -> 'FloatArray':
        """ Tidal susceptibility [N m] """
        return self._tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def unique_tidal_frequencies(self) -> Dict['FreqSig', 'FloatArray']:
        """ Unique tidal frequencies (abs(tidal modes)) stored by frequency signature """
        return self._unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def tidal_terms_by_frequency(self) -> Dict['FreqSig', Dict[int, 'DissipTermsArray']]:
        """ Tidal terms stored by frequency signature """
        return self._tidal_terms_by_frequency

    @tidal_terms_by_frequency.setter
    def tidal_terms_by_frequency(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def tidal_heating_global(self) -> 'FloatArray':
        """ Global tidal heating rate [W] """
        return self._tidal_heating_global

    @tidal_heating_global.setter
    def tidal_heating_global(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def global_negative_imk_by_orderl(self) -> Dict[int, 'FloatArray']:
        """ Global negative of the imaginary portion of the Love number, -Im[k_l] """
        return self._global_negative_imk_by_orderl

    @global_negative_imk_by_orderl.setter
    def global_negative_imk_by_orderl(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def global_love_by_orderl(self) -> Dict[int, 'FloatArray']:
        """ Global complex Love number, k_l """
        return self._global_love_by_orderl

    @global_love_by_orderl.setter
    def global_love_by_orderl(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def effective_q_by_orderl(self) -> Dict[int, 'FloatArray']:
        """ World's effective tidal dissipation factor (for each tidal order level) """
        return self._effective_q_by_orderl

    @effective_q_by_orderl.setter
    def effective_q_by_orderl(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def dUdM(self) -> 'FloatArray':
        """ Global partial derivative of the tidal potential with respect to the mean anomaly """
        return self._dUdM

    @dUdM.setter
    def dUdM(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def dUdw(self) -> 'FloatArray':
        """ Global partial derivative of the tidal potential with respect to the argument of pericentre """
        return self._dUdw

    @dUdw.setter
    def dUdw(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def dUdO(self) -> 'FloatArray':
        """ Global partial derivative of the tidal potential with respect to the argument of node """
        return self._dUdO

    @dUdO.setter
    def dUdO(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def fixed_q(self) -> float:
        """ Fixed tidal dissipation factor (used in CPL tidal calculation method) """
        return self._fixed_q

    @fixed_q.setter
    def fixed_q(self, new_fixed_q: float):
        self.set_fixed_q(new_fixed_q)

    @property
    def fixed_dt(self) -> float:
        """ Fixed tidal dissipation frequency dependency (used in CTL tidal calculation method) """
        return self._fixed_dt

    @fixed_dt.setter
    def fixed_dt(self, new_fixed_dt: float):
        self.set_fixed_dt(new_fixed_dt)

    @property
    def fixed_k2(self) -> float:
        """ Fixed static k2 for the world """
        return self._fixed_k2

    @fixed_k2.setter
    def fixed_k2(self, value):
        raise OuterscopePropertySetError

    # # Outer-scope Properties
    # World properties
    @property
    def semi_major_axis(self):
        """ Outer-scope wrapper for world.semi_major_axis """
        return self.world.semi_major_axis

    @semi_major_axis.setter
    def semi_major_axis(self, value):
        raise OuterscopePropertySetError

    @property
    def orbital_frequency(self):
        """ Outer-scope wrapper for world.orbital_frequency """
        return self.world.orbital_frequency

    @orbital_frequency.setter
    def orbital_frequency(self, value):
        raise OuterscopePropertySetError

    @property
    def spin_frequency(self):
        """ Outer-scope wrapper for world.spin_frequency """
        return self.world.spin_frequency

    @spin_frequency.setter
    def spin_frequency(self, value):
        raise OuterscopePropertySetError

    @property
    def eccentricity(self):
        """ Outer-scope wrapper for world.eccentricity """
        return self.world.eccentricity

    @eccentricity.setter
    def eccentricity(self, value):
        raise OuterscopePropertySetError

    @property
    def obliquity(self):
        """ Outer-scope wrapper for world.obliquity """
        return self.world.obliquity

    @obliquity.setter
    def obliquity(self, value):
        raise OuterscopePropertySetError

    @property
    def tidal_host(self):
        """ Outer-scope wrapper for world.orbit.tidal_host """
        return self.world.tidal_host

    @tidal_host.setter
    def tidal_host(self, value):
        raise OuterscopePropertySetError

    @property
    def radius(self):
        """ Outer-scope wrapper for world.radius """
        return self.world.radius

    @radius.setter
    def radius(self, value):
        raise OuterscopePropertySetError

    @property
    def moi(self):
        """ Outer-scope wrapper for world.moi """
        return self.world.moi

    @moi.setter
    def moi(self, value):
        raise OuterscopePropertySetError

    # # Dunder methods
    def __str__(self):
        str_ = f'{self.__class__.__name__}'
        if self.world is not None:
            str_ += f' ({self.world})'

        return str_

    # # Static Methods
    @staticmethod
    def calculate_tidal_susceptibility(
        host_mass: float, target_radius: float,
        semi_major_axis: 'FloatArray'
        ) -> 'FloatArray':
        """ Calculate the tidal susceptibility for a target object orbiting

        Wrapper for TidalPy.tides.dissipation.calc_tidal_susceptibility

        Parameters
        ----------
        host_mass : float
            Tidal host's mass [kg]
        target_radius : float
            Target body's mean radius [m]
        semi_major_axis : FloatArray
            Orbital separation between the target and host [m]

        Returns
        -------
        tidal_susceptibility : FloatArray
            Tidal Susceptibility [N m]
        """

        tidal_susceptibility = calc_tidal_susceptibility(host_mass, target_radius, semi_major_axis)
        return tidal_susceptibility

    @staticmethod
    def calculate_effective_rigidity(
        shear_modulus: 'FloatArray', gravity: float, radius: float, bulk_density: float,
        tidal_order_l: int = 2
        ) -> 'FloatArray':
        """ Calculate the effective rigidity of a layer or planet

        Wrapper for TidalPy.tides.love1d.effective_rigidity_general

        Parameters
        ----------
        shear_modulus : FloatArray
            Shear modulus of a layer/planet [Pa]
        gravity : float
            Acceleration of gravity at the surface of a layer/planet [m s-2]
        radius : float
            Radius at the top of a layer/planet [m]
        bulk_density : float
            Bulk density of a layer/planet [kg m-3]
        tidal_order_l : int = 2
            Tidal harmonic order integer

        Returns
        -------
        effective_rigidity : FloatArray
            Effective rigidity of the layer/planet [Pa Pa-1]
        """

        effective_rigidity = effective_rigidity_general(
            shear_modulus, gravity, radius, bulk_density,
            order_l=tidal_order_l
            )
        return effective_rigidity

    @staticmethod
    def calculate_complex_love_number(
        shear_modulus: 'FloatArray', complex_compliance: 'ComplexArray',
        effective_rigidity: 'FloatArray', tidal_order_l: int = 2
        ) -> 'ComplexArray':
        """ Calculate the complex Love number of a layer or planet

        Wrapper for TidalPy.tides.love1d.complex_love_general

        Parameters
        ----------
        shear_modulus : FloatArray
            Shear modulus of a layer/planet [Pa]
        complex_compliance : ComplexArray
            Complex compliance of a layer/planet [Pa -1]
        effective_rigidity : FloatArray
            Effective rigidity of a layer/planet [Pa Pa-1]
            See: Tides.calculate_effective_rigidity
        tidal_order_l : int = 2
            Tidal harmonic order integer

        Returns
        -------
        complex_love_number : ComplexArray
            Complex Love number for the layer/planet
        """

        complex_love_number = complex_love_general(
            complex_compliance, shear_modulus, effective_rigidity,
            order_l=tidal_order_l
            )
        return complex_love_number
