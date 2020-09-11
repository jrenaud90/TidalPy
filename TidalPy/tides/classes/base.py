""" Base Tides Module
"""

from typing import TYPE_CHECKING, Dict, Tuple

import numpy as np

from exceptions import (AttributeNotSetError, ImproperPropertyHandling, OuterscopePropertySetError,
                        ConfigPropertyChangeError, FailedForcedStateUpdate,
                        ImplementedBySubclassError, IncompatibleModelError)
from .defaults import tide_defaults
from ..dissipation import calc_tidal_susceptibility, calc_tidal_susceptibility_reduced
from ..eccentricityFuncs import EccenOutput
from ..inclinationFuncs import InclinOutput
from ..love1d import complex_love_general, effective_rigidity_general
from ..mode_manipulation import find_mode_manipulators, FreqSig, DissipTermsArray
from ... import log
from ...utilities.classes.config.config import WorldConfigHolder
from ...utilities.types import FloatArray, ComplexArray

if TYPE_CHECKING:
    from ...structures.worlds import TidalWorldType

# TODO: Add a spin-sync version


class TidesBase(WorldConfigHolder):

    """ TidesBase Class - Holder for all tidal heating and tidal potential calculations

    Tides class stores model parameters and methods for calculating tidal heating and tidal potential derivatives
        which are general functions of (T, P, melt_frac, w, e, obliquity)

    Attributes
    ----------
    _thermal_set : bool
        Flag - If world's thermal state (can calculate complex compliance) has been set or not.
        Class Property
    _orbit_set : bool
        Flag - If the world's orbital state has been set or not.
        Class Property
    _tidal_susceptibility : np.ndarray
        Tidal susceptibility for the world [N m].
        Class Property
    _tidal_susceptibility_reduced : np.ndarray
        Reduced tidal susceptibility for the world (no semi-major axis dependence).
        Class Property
    _unique_tidal_frequencies : Dict[freq_sig, np.ndarray]
        Unique tidal forcing frequencies experienced by the world
            (combination of orbital motion and spin frequency).
        Class Property
    _tidal_terms_by_frequency : Dict[freq_sig, Dict[int, np.ndarray]]

        Class Property
    """

    model = 'base'
    default_config = tide_defaults['base']
    world_config_key = 'tides'

    def __init__(self, world: 'TidalWorldType', store_config_in_world: bool = True, auto_compile_funcs: bool = True):

        super().__init__(world, store_config_in_world=store_config_in_world)

        # Flags
        self._thermal_set = False
        self._orbit_set = False

        # State properties
        self._tidal_susceptibility = None
        self._tidal_susceptibility_reduced = None
        self._unique_tidal_frequencies = None
        self._tidal_terms_by_frequency = None
        self._tidal_heating_global = None
        self._negative_imk_global = None
        self._dUdM = None
        self._dUdw = None
        self._dUdO = None
        self._tidal_polar_torque = None
        self._spin_rate_derivative = None

        # Model configurations that will be set in setup
        self._eccentricity_truncation_lvl = None
        self._max_tidal_order_lvl = None
        self._use_obliquity_tides = None

        # Functions to be initialized in self.setup
        self._eccentricity_results = None
        self._obliquity_results = None
        self.eccentricity_func = None
        self.obliquity_func = None
        self.collapse_modes_func = None
        self.calculate_modes_func = None

        # TidalPy logging and debug info
        log.debug(f'Building {self.model} tides class for {self.world.name}.')

        # Call setup for initialization
        self.setup(auto_compile_funcs=auto_compile_funcs)

    def setup(self, overload_tidal_l: int = None, overload_eccentricity_truncation: int = None,
              auto_compile_funcs: bool = False):
        """ Load configurations into the Tides class and import any config-dependent functions.

        This setup process is separate from the __init__ method because the Orbit class may need to overload some
            configurations after class initialization.
        """

        # Load in configurations
        self._use_obliquity_tides = self.config['obliquity_tides_on']
        if overload_tidal_l is not None:
            log(f'Tidal-l overloaded in {self}', level='debug')
            self._max_tidal_order_lvl = overload_tidal_l
        else:
            self._max_tidal_order_lvl = self.config['max_tidal_order_l']

        if overload_eccentricity_truncation is not None:
            log(f'Eccentricity Truncation overloaded in {self}', level='debug')
            self._eccentricity_truncation_lvl = overload_eccentricity_truncation
        else:
            self._eccentricity_truncation_lvl = self.config['eccentricity_truncation_lvl']

        # Setup functions
        self.eccentricity_func, self.obliquity_func, self.collapse_modes_func, self.calculate_modes_func = \
            find_mode_manipulators(self.max_tidal_order_lvl, self.eccentricity_truncation_lvl, self.use_obliquity_tides)

        # Clear the state as some of the old values probably would change with an override update
        self.clear_state()

    def clear_state(self):

        super().clear_state()

        self._eccentricity_results = None
        self._obliquity_results = None
        self._tidal_susceptibility = None
        self._tidal_susceptibility_reduced = None
        self._unique_tidal_frequencies = None
        self._tidal_terms_by_frequency = None
        self._tidal_heating_global = None
        self._negative_imk_global = None
        self._dUdM = None
        self._dUdw = None
        self._dUdO = None
        self._tidal_polar_torque = None
        self._spin_rate_derivative = None

    def initialize_tides(self):
        """ Initialize various tidal parameters once a tidal host is connected to the target body.
        """

        if self.tidal_host is None:
            raise AttributeNotSetError('Tidal host must be connected to target body in order to initialize tides.')

        # Try to update properties that depend on orbit and tidal host
        self._tidal_susceptibility_reduced = \
            calc_tidal_susceptibility_reduced(self.tidal_host.mass, self.world.radius)

        if self.semi_major_axis is not None:
            self._tidal_susceptibility = calc_tidal_susceptibility(self.tidal_host.mass, self.world.radius,
                                                                   self.semi_major_axis)

    def update_orbit_spin(self, force_update: bool = True, eccentricity_change: bool = True,
                          obliquity_change: bool = True, frequency_change: bool = True,
                          eccentricity_results: Dict[int, EccenOutput] = None) -> \
            Tuple[Dict[FreqSig, np.ndarray], Dict[FreqSig, Dict[int, DissipTermsArray]]]:
        """ Calculate tidal heating and potential derivative terms based on the current orbital state.

        This will also calculate new unique tidal frequencies which must then be digested by the rheological model
            at each planetary layer.


        Parameters
        ----------
        force_update : bool = True
            If True, will raise an error if the update does not successful complete.
            Failed completions are usually a result of a missing state property(ies).
        eccentricity_change : bool = True
            If there was no change in eccentricity (or if the orbit set the eccentricity) set this to False for a
            performance boost. If False, eccentricity functions won't be called.
        obliquity_change : bool = True
            If there was no change in obliquity set this to False for a performance boost.
            If False, obliquity functions won't be called.
        frequency_change : bool = True
            If there was no change in orbital or rotation frequency set this to False for a performance boost.
            If False, calculate_terms won't be called.
        eccentricity_results : Dict[int, EccenOutput] = None
            Eccentricity functions can be calculated by the Orbit class (or by the user). Pass a pre-caclualted result
            for a performance boost.

        Returns
        -------
        unique_tidal_frequencies : Dict[FreqSig, np.ndarray]
            Each unique frequency stored as a signature (orbital motion and spin-rate coeffs), and the calculated frequency
                (combo of orbital motion and spin-rate) [rad s-1]
        tidal_terms_by_frequency : Dict[FreqSig, Dict[int, DissipTermsArray]]
            Results for tidal heating, dU/dM, dU/dw, dU/dO are stored in a tuple for each tidal harmonic l and
                unique frequency.

        See Also
        --------
        TidalPy.tides.dissipation.mode_collapse
        """

        # Check to see if all the needed state properties are present
        all_values_present = True
        for param in [self.orbital_frequency, self.spin_frequency, self.eccentricity, self.semi_major_axis]:
            if param is None:
                all_values_present = False
                break

        if all_values_present:
            need_to_collapse_modes = False

            if eccentricity_results is not None:
                self._eccentricity_results = eccentricity_results
                need_to_collapse_modes = True
            elif eccentricity_change:
                self._eccentricity_results = self.eccentricity_func(self.eccentricity)
                need_to_collapse_modes = True

            if obliquity_change:
                self._obliquity_results = self.obliquity_func(self.obliquity)
                need_to_collapse_modes = True

            # Check that obliquity and eccentricity results have the same length (same max_l used).
            if len(self.obliquity_results) != len(self.eccentricity_results):
                raise IncompatibleModelError('Obliquity and Eccentricity results do not have the same length.'
                                             'max_tidal_l may not have been set the same for the functions.')

            if frequency_change:
                self._unique_tidal_frequencies, self._tidal_terms_by_frequency = \
                    self.calculate_modes_func(self.spin_frequency, self.orbital_frequency, self.semi_major_axis,
                                              self.radius, self.eccentricity_results, self.obliquity_results)
                # Now that there are new frequencies, tell the world so that new complex compliances can be calculated.
                self.world.unique_frequencies_updated()
                need_to_collapse_modes = True

            if frequency_change:
                # Orbital changes may have changed the tidal susceptibility
                self._tidal_susceptibility = calc_tidal_susceptibility(self.tidal_host.mass, self.world.radius,
                                                                       self.semi_major_axis)

            # Set orbit to true. Check if thermals are set and collapse modes
            self._orbit_set = True
            if self.thermal_set:
                self.collapse_modes(force_update=force_update)

        else:
            if force_update:
                raise FailedForcedStateUpdate

        # Return frequencies and tidal terms
        return self.unique_tidal_frequencies, self.tidal_terms_by_frequency

    def collapse_modes(self, force_update: bool = True) -> DissipTermsArray:
        """ Calculate Global Love number based on current thermal state.

        Requires a prior update_orbit_spin() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        Parameters
        ----------
        force_update : bool = True
            If True, will raise an error if the update does not successful complete.
            Failed completions are usually a result of a missing state property(ies).

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
        TidalPy.tides.Tides.update_orbit_spin
        """

        raise ImplementedBySubclassError('mode_collapse is implemented by a child class of TidesBase')

    def set_spin_derivative(self) -> np.ndarray:
        """ Calculate spin-rate derivative based on current state

        Requires a prior thermal_change() call as dUdO must be set before spin-rate derivative can be calculated
            is called.

        Returns
        -------
        spin_rate_derivative : np.ndarray
            Spin-rate derivative for the planet [rads s-2]
        """

        if self.dUdO is None:
            raise AttributeNotSetError(f'Potential derivatives not calculated for {self.world.name}.')

        spin_rate_derivative = self.tidal_host.mass * self.dUdO / self.moi
        self._spin_rate_derivative = spin_rate_derivative
        self._tidal_polar_torque = self.tidal_host.mass * self.dUdO

        return spin_rate_derivative

    @staticmethod
    def calculate_tidal_susceptibility(host_mass: float, target_radius: float, semi_major_axis: FloatArray) \
            -> FloatArray:
        """ Calculate the tidal susceptibility for a target object orbiting

        Wrapper for dissipation.py/calc_tidal_susceptibility

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
    def calculate_effective_rigidity(shear_modulus: FloatArray, gravity: float, radius: float, bulk_density: float,
                                     tidal_order_l: int = 2) -> FloatArray:
        """ Calculate the effective rigidity of a layer or planet

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

        effective_rigidity = effective_rigidity_general(shear_modulus, gravity, radius, bulk_density,
                                                        order_l=tidal_order_l)
        return effective_rigidity

    @staticmethod
    def calculate_complex_love_number(shear_modulus: FloatArray, complex_compliance: ComplexArray,
                                      effective_rigidity: FloatArray, tidal_order_l: int = 2) -> ComplexArray:
        """ Calculate the complex Love number of a layer or planet

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

        complex_love_number = complex_love_general(complex_compliance, shear_modulus, effective_rigidity,
                                                   order_l=tidal_order_l)
        return complex_love_number


    # Configuration properties
    @property
    def eccentricity_truncation_lvl(self) -> int:
        return self._eccentricity_truncation_lvl

    @eccentricity_truncation_lvl.setter
    def eccentricity_truncation_lvl(self, value):
        # TODO: Think about if you want the user to update these. These setters could make a call to self.setup()
        #    which the user could make on their own. So it may make sense to allow the setter.
        raise ConfigPropertyChangeError

    @property
    def max_tidal_order_lvl(self) -> int:
        return self._max_tidal_order_lvl

    @max_tidal_order_lvl.setter
    def max_tidal_order_lvl(self, value):
        raise ConfigPropertyChangeError

    @property
    def use_obliquity_tides(self) -> bool:
        return self._use_obliquity_tides

    @use_obliquity_tides.setter
    def use_obliquity_tides(self, value):
        raise ConfigPropertyChangeError


    # State properties
    @property
    def thermal_set(self) -> bool:
        return self._thermal_set

    @thermal_set.setter
    def thermal_set(self, value):
        raise ImproperPropertyHandling

    @property
    def orbit_set(self) -> bool:
        return self._orbit_set

    @orbit_set.setter
    def orbit_set(self, value):
        raise ImproperPropertyHandling

    @property
    def eccentricity_results(self) -> Dict[int, InclinOutput]:
        return self._eccentricity_results

    @eccentricity_results.setter
    def eccentricity_results(self, value):
        raise ImproperPropertyHandling

    @property
    def obliquity_results(self) -> Dict[int, EccenOutput]:
        return self._obliquity_results

    @obliquity_results.setter
    def obliquity_results(self, value):
        raise ImproperPropertyHandling

    @property
    def tidal_susceptibility_reduced(self) -> np.ndarray:
        return self._tidal_susceptibility_reduced

    @tidal_susceptibility_reduced.setter
    def tidal_susceptibility_reduced(self, value):
        raise ImproperPropertyHandling

    @property
    def tidal_susceptibility(self) -> np.ndarray:
        return self._tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise ImproperPropertyHandling

    @property
    def unique_tidal_frequencies(self) -> Dict[FreqSig, np.ndarray]:
        return self._unique_tidal_frequencies

    @unique_tidal_frequencies.setter
    def unique_tidal_frequencies(self, value):
        raise ImproperPropertyHandling

    @property
    def tidal_terms_by_frequency(self) -> Dict[FreqSig, Dict[int, DissipTermsArray]]:
        return self._tidal_terms_by_frequency

    @tidal_terms_by_frequency.setter
    def tidal_terms_by_frequency(self, value):
        raise ImproperPropertyHandling

    @property
    def tidal_heating_global(self) -> np.ndarray:
        return self._tidal_heating_global

    @tidal_heating_global.setter
    def tidal_heating_global(self, value):
        raise ImproperPropertyHandling

    @property
    def negative_imk_global(self) -> np.ndarray:
        return self._negative_imk_global

    @negative_imk_global.setter
    def negative_imk_global(self, value):
        raise ImproperPropertyHandling

    @property
    def dUdM(self) -> np.ndarray:
        return self._dUdM

    @dUdM.setter
    def dUdM(self, value):
        raise ImproperPropertyHandling

    @property
    def dUdw(self) -> np.ndarray:
        return self._dUdw

    @dUdw.setter
    def dUdw(self, value):
        raise ImproperPropertyHandling

    @property
    def dUdO(self) -> np.ndarray:
        return self._dUdO

    @dUdO.setter
    def dUdO(self, value):
        raise ImproperPropertyHandling

    @property
    def tidal_polar_torque(self) -> np.ndarray:
        return self._tidal_polar_torque

    @tidal_polar_torque.setter
    def tidal_polar_torque(self, value):
        raise ImproperPropertyHandling

    @property
    def spin_rate_derivative(self) -> np.ndarray:
        return self._spin_rate_derivative

    @spin_rate_derivative.setter
    def spin_rate_derivative(self, value):
        raise ImproperPropertyHandling


    # Outer-scope Properties
    @property
    def semi_major_axis(self):
        return self.world.semi_major_axis

    @semi_major_axis.setter
    def semi_major_axis(self, value):
        raise OuterscopePropertySetError

    @property
    def orbital_frequency(self):
        return self.world.orbital_frequency

    @orbital_frequency.setter
    def orbital_frequency(self, value):
        raise OuterscopePropertySetError

    @property
    def spin_frequency(self):
        return self.world.spin_frequency

    @spin_frequency.setter
    def spin_frequency(self, value):
        raise OuterscopePropertySetError

    @property
    def eccentricity(self):
        return self.world.eccentricity

    @eccentricity.setter
    def eccentricity(self, value):
        raise OuterscopePropertySetError

    @property
    def obliquity(self):
        if self.use_obliquity_tides:
            return self.world.obliquity
        else:
            return np.zeros_like(self.eccentricity)

    @obliquity.setter
    def obliquity(self, value):
        raise OuterscopePropertySetError

    @property
    def tidal_host(self):
        return self.world.tidal_host

    @tidal_host.setter
    def tidal_host(self, value):
        raise OuterscopePropertySetError

    @property
    def radius(self):
        return self.world.radius

    @radius.setter
    def radius(self, value):
        raise OuterscopePropertySetError

    @property
    def moi(self):
        return self.world.moi

    @moi.setter
    def moi(self, value):
        raise OuterscopePropertySetError

