from __future__ import annotations

from typing import List, Union

import numpy as np
from scipy.constants import G

from TidalPy.utilities.numpyHelper.array_other import value_cleanup
from .. import debug_mode
from ..dynamics import diff_eqs_duel_dissipation, diff_eqs_single_dissipation
from ..exceptions import (ImproperPropertyHandling, ArgumentException, IncompatibleModelConfigError,
                          IncorrectArgumentType, IncorrectModelInitialized, ParameterValueError, ParameterMissingError)
from ..initialize import log
from ..structures.worlds import BasicWorld, Star, TidalWorld, WorldBase
from TidalPy.utilities.types import FloatArray
from ..utilities.classes import TidalPyClass
from TidalPy.tools.conversions import Au2m, days2rads, rads2days


HostTypes = Union[BasicWorld, TidalWorld]
PlanetRefType = Union[str, int, WorldBase]
TargetBodyType = Union[TidalWorld, List[TidalWorld]]


def pull_out_orbit_defaults(planet_obj):

    # Update those dummy variables if information was provided in the configurations
    orbital_freq = planet_obj.config.get('orbital_freq', None)
    orbital_period = planet_obj.config.get('orbital_period', None)
    semi_major_axis = planet_obj.config.get('semi_major_axis', None)
    semi_major_axis_inau = planet_obj.config.get('semi_major_axis_in_au', False)
    eccentricity = planet_obj.config.get('eccentricity', None)
    inclination = planet_obj.config.get('inclination', None)

    if orbital_period is not None:
        if orbital_freq is not None:
            log(f'Both orbital frequency and period were provided for {planet_obj.name}. '
                f'Using frequency instead.', level='info')
        else:
            orbital_freq = np.asarray([2. * np.pi / (orbital_period * 24. * 60. * 60.)])
    if semi_major_axis is not None:
        semi_major_axis = np.asarray([semi_major_axis])
        if semi_major_axis_inau:
            semi_major_axis = Au2m(semi_major_axis)
    if eccentricity is not None:
        eccentricity = np.asarray([eccentricity])
    if inclination is not None:
        inclination = np.asarray([inclination])

    return orbital_freq, semi_major_axis, eccentricity, inclination


class OrbitBase(TidalPyClass):
    """ OrbitBase class connects stars, hosts, and target planets---handling all mutual calculations

    """

    def __init__(self, star: Star, host: HostTypes = None, target_bodies: TargetBodyType = None,
                 duel_dissipation: bool = False, host_tide_raiser_location: int = None, use_host_orbit: bool = False,
                 time_study: bool = False):
        """

        Parameters
        ----------
        star : Star
        host : HostTypes
        target_bodies : TargetBodyType
        duel_dissipation : bool
        host_tide_raiser_location : int
        time_study : bool
        """

        super().__init__()

        # Load in flags
        self.duel_dissipation = duel_dissipation
        self.is_time_study = time_study

        # Load in planets and stars
        self.star = star
        if host is None:
            host = star
        self.host = host
        self.host.is_host = True
        if target_bodies is None:
            self.target_bodies = list()
        else:
            if isinstance(target_bodies, TidalWorld):
                self.target_bodies = [target_bodies]  # type: List[TidalWorld]
            elif type(target_bodies) == list:
                self.target_bodies = target_bodies    # type: List[TidalWorld]
            else:
                raise TypeError
        self.all_objects = [self.star, self.host] + self.target_bodies
        self.all_objects_byname = {self.star.name: self.star}
        if self.host.name not in self.all_objects_byname:
            self.all_objects_byname[self.host.name] = self.host
        for target_body in self.target_bodies:
            self.all_objects_byname[target_body.name] = target_body

        # Load in parameters into planets
        for orbital_loc, world in enumerate(self):
            world._orbit = self
            world.orbit_location = orbital_loc

        # Store aliases
        need_to_add = []
        for name, obj in self.all_objects_byname.items():
            if name.lower() not in self.all_objects_byname:
                need_to_add.append((name.lower(), obj))
            if name.title() not in self.all_objects_byname:
                need_to_add.append((name.title(), obj))
        for new_name, obj in need_to_add:
            self.all_objects_byname[new_name] = obj

        # Perform some checks
        if self.duel_dissipation and not isinstance(self.host, TidalWorld):
            raise IncompatibleModelConfigError('Duel dissipation requires a tidal (TidalWorld class) host.')

        # The host body may be tidally dissipating due to one of the other target planets.
        self._host_tide_raiser_loc = None
        if host_tide_raiser_location is not None:
            self.host_tide_raiser_loc = host_tide_raiser_location
        elif len(self.target_bodies) == 1:
            self.host_tide_raiser_loc = 2
        elif len(self.target_bodies) > 1:
            log.warn('More than one target body provided. Assuming that innermost is the primary tide raiser on host.')
            self.host_tide_raiser_loc = 2

        # State orbital variables (must be at least 2: for the star and the host)
        self._eccentricities = [np.asarray(0.), np.asarray(0.)]  # type: List[Union[None, np.ndarray]]
        self._inclinations = [np.asarray(0.), np.asarray(0.)]    # type: List[Union[None, np.ndarray]]
        self._orbital_freqs = [None, None]                       # type: List[Union[None, np.ndarray]]
        self._semi_major_axis = [None, None]                     # type: List[Union[None, np.ndarray]]
        self._derivative_a = [np.asarray(0.), np.asarray(0.)]    # type: List[Union[None, np.ndarray]]
        self._derivative_e = [np.asarray(0.), np.asarray(0.)]    # type: List[Union[None, np.ndarray]]

        # Are target bodies orbiting the star or the host
        self.star_host = False
        if self.star is self.host:
            self.star_host = True
        self.equilib_distance_funcs = [None, lambda: self.star.semi_major_axis]
        self.equilib_eccentricity_funcs = [None, lambda: self.star.eccentricity]

        # Is the host a tidal body?
        self.has_tidal_host = isinstance(self.host, TidalWorld)

        # Determine various constants that depend upon other objects properties
        for target_body in self.target_bodies:

            # Check for potential issues
            if target_body is self.star:
                raise ParameterValueError("The orbit's star can not be a target body")
            if target_body is self.host:
                raise ParameterValueError("The orbit's host can not be a target body. "
                                     "Nevertheless, tides can still be calculated in the host. See documentation.")

            # Equilibrium Temperature
            if self.star_host:
                self.equilib_distance_funcs.append(lambda: target_body.semi_major_axis)
                self.equilib_eccentricity_funcs.append(lambda: target_body.eccentricity)
            else:
                self.equilib_distance_funcs.append(lambda: self.star.semi_major_axis)
                self.equilib_eccentricity_funcs.append(lambda: self.star.eccentricity)

            # Inflated Tidal Susceptibility
            target_body.tidal_susceptibility_inflated = (3. / 2.) * G * self.host.mass**2 * target_body.radius**5

            # Store dummy values for the state variables
            self._eccentricities.append(np.asarray([0.]))
            self._inclinations.append(np.asarray([0.]))
            self._orbital_freqs.append(None)
            self._semi_major_axis.append(None)
            self._derivative_a.append(np.asarray([0.]))
            self._derivative_e.append(np.asarray([0.]))

        for t_i, world in enumerate(self.all_objects):
            # Star does not need (or have) any of the following parameters - so skip it
            if world is self.star:
                continue

            # Update those dummy variables if information was provided in the configurations
            orbital_freq, semi_major_axis, eccentricity, inclination = pull_out_orbit_defaults(world)

            if not self.star_host and world is self.host:
                # Initial parameters set in the config file are assumed to be set from star position
                self.set_orbit(self.star, orbital_freq, semi_major_axis, eccentricity, inclination,
                               set_by_planet=False, force_calculation=False)

            if not self.star_host and use_host_orbit and world is not self.host:
                # Use host's orbital parameters instead
                orbital_freq_host, semi_major_axis_host, eccentricity_host, inclination_host = \
                    pull_out_orbit_defaults(self.host)
                self.set_orbit(world, orbital_freq_host, semi_major_axis_host, eccentricity_host, inclination_host,
                               set_by_planet=False, force_calculation=False)
            else:
                self.set_orbit(world, orbital_freq, semi_major_axis, eccentricity, inclination,
                               set_by_planet=False, force_calculation=False)

        # Attempt to initialize insolation heating
        for target_body in self.target_bodies:
            self.calculate_insolation(target_body)
        if isinstance(self.host, TidalWorld):
            self.calculate_insolation(self.host)

    def clear_state(self):

        if debug_mode:
            log(f'State cleared for {self}', level='debug')

        # Clear state properties
        # TODO:

        # Clear planet properties
        for world in self:
            # We are 'preserving' the orbit here because all of those properties should have just been reset by the
            #    above section.
            world.clear_state(preserve_orbit=True)

    def find_planet_pointer(self, planet_reference: PlanetRefType) -> WorldBase:
        """ Find the object pointer to a planet or object stored in this orbit

        Possible references:
            str:
                Name of the planet provided in the configuration
                    This is case insensitive unless there the config provided a capital letter any where other than
                    the first letter.
            int:
                Planet's location in the orbit based on the following scheme:
                    0 == Star
                    1 == Host
                    2+ == Various target planets (these are probably the ones you want to set!)

        Parameters
        ----------
        planet_reference : PlanetRefType
            User-friendly reference to a planet or object stored in this orbit

        Returns
        -------
        planet_reference : TidalWorld
            Pointer to the planet object
        """

        if isinstance(planet_reference, WorldBase):
            # Already is the pointer! Make sure it is actually loaded into this orbit.

            for obj in self.all_objects:
                if planet_reference is obj:
                    break
            else:
                # pointer is not actually a pointer for any object that this orbit was initialized for.
                raise ArgumentException('Planet Reference Pointer not found in orbit instance.')

            return planet_reference

        if type(planet_reference) == int:
            return self.all_objects[planet_reference]
        if type(planet_reference) == str:
            try:
                return self.all_objects_byname[planet_reference]
            except KeyError:
                return self.all_objects_byname[planet_reference.lower()]
        raise IncorrectArgumentType

    # Functionality to set the orbit state variables
    def set_orbit(self, planet_reference: PlanetRefType, orbital_freq: FloatArray = None,
                  semi_major_axis: FloatArray = None, eccentricity: FloatArray = None, inclination: FloatArray = None,
                  semi_major_axis_in_au: bool = False, inclination_in_deg: bool = False,
                  set_by_planet: bool = False, force_calculation: bool = True):
        """ Set the orbital state (orbital frequency, eccentricity, etc.) of a planet in this orbit.

        Can set orbital frequency, semi-major axis, eccentricity, and/or inclination of a planet at the
        orbit location: planet_loc. Use this method instead of the individual set_eccentricity, set_inclination, etc.
        methods when you are changing 2 or more parameters simultaneously (it is more efficient).

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference used to find the planet
        orbital_freq : FloatArray
            Orbital frequency in [rads s-1]
            Optional, Mutually exclusive with semi-major axis
        semi_major_axis : FloatArray
            Semi-major axis in [m]
            Optional, Mutually exclusive with orbital_freq
        eccentricity : FloatArray
            Eccentricity (only elliptical and circular orbits are supported in TidalPy)
            Optional
        inclination : FloatArray
            Orbital inclination relative to the orbital plane in [rads]
            Optional
        semi_major_axis_in_au : bool
            If True then TidalPy will convert semi-major axis to meters (unit required by calculations)
        inclination_in_deg : bool = False
            If true then TidalPy will convert the input into radians (the unit that is required for calculations)
        set_by_planet : bool = False
            A flag used by TidalPy planets to make calls to this method. Leave it False for user use.
            If set to true then other models and parameters will not update correctly.
        force_calculation : bool = True
            If set to true then ParameterMissingErrors will be raised

        See Also
        --------
        OrbitBase.set_orbital_freq
        OrbitBase.set_semi_major_axis
        OrbitBase.set_eccentricity
        OrbitBase.set_inclination
        """

        planet_pointer = self.find_planet_pointer(planet_reference)

        # The "set_by_planet" are set to True below even if the "set_by_planet" is False above because the function
        #    Will do a single call to orbit_update at the end.
        if orbital_freq is not None:
            if semi_major_axis is not None:
                raise ImproperPropertyHandling('Only set orbital frequency or semi-major axis, not both.')
            self.set_orbital_freq(planet_pointer, orbital_freq, set_by_planet=True)

        if semi_major_axis is not None:
            self.set_semi_major_axis(planet_pointer, semi_major_axis, semi_major_axis_in_au=semi_major_axis_in_au,
                                     set_by_planet=True)

        if eccentricity is not None:
            self.set_eccentricity(planet_pointer, eccentricity, set_by_planet=True)

        if inclination is not None:
            self.set_inclination(planet_pointer, inclination, inclination_in_deg=inclination_in_deg, set_by_planet=True)

        if planet_pointer is not self.star:
            try:
                self.update_orbit(planet_reference, set_by_planet=set_by_planet)
            except ParameterMissingError as error:
                if force_calculation:
                    raise error

    def set_orbital_freq(self, planet_reference: PlanetRefType, new_orbital_freq: FloatArray,
                         set_by_planet: bool = False):
        """ Set the orbital frequency of a planet at planet_loc

        Use Orbit.set_orbit if setting more than one state parameter.

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference used to find the planet
        new_orbital_freq : FloatArray
            Orbital frequency in [rads s-1]
            Optional, Mutually exclusive with semi-major axis
        set_by_planet : bool = False
            A flag used by TidalPy planets to make calls to this method. Leave it False for user use.
            If set to true then other models and parameters will not update correctly.

        See Also
        --------
        OrbitBase.set_orbit
        """

        planet_pointer = self.find_planet_pointer(planet_reference)
        planet_loc = planet_pointer.orbit_location

        if planet_pointer is self.host and not self.star_host and not set_by_planet:
            if debug_mode:
                # This would change the orbit between the host and the star, not between the target body and the star.
                # This may not be what the user wants to do.
                log("Attempting to change the host planet's orbit", level='debug')

        self._orbital_freqs[planet_loc] = value_cleanup(new_orbital_freq)

        # Changing the orbital frequency also changes the semi-major axis. Update via Kepler Laws
        if planet_pointer is self.host:
            self._semi_major_axis[planet_loc] = \
                np.cbrt(G * (planet_pointer.mass + self.star.mass) / new_orbital_freq**2)
        else:
            # Star and Target bodies all orbit the host (albeit the star doesn't do much...).
            self._semi_major_axis[planet_loc] = \
                np.cbrt(G * (planet_pointer.mass + self.host.mass) / new_orbital_freq**2)
        if not set_by_planet:
            self.update_orbit(planet_reference, set_by_planet=False)

    def set_semi_major_axis(self, planet_reference: PlanetRefType, new_semi_major_axis: FloatArray,
                            semi_major_axis_in_au: bool = False, set_by_planet: bool = False):
        """ Set the semi-major axis of a planet at planet_loc

        Use Orbit.set_orbit if setting more than one state parameter.

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference used to find the planet
        new_semi_major_axis : FloatArray
            Semi-major axis in [m]
            Optional, Mutually exclusive with orbital_freq
        semi_major_axis_in_au : bool
            If True then TidalPy will convert semi-major axis to meters (unit required by calculations)
        set_by_planet : bool = False
            A flag used by TidalPy planets to make calls to this method. Leave it False for user use.
            If set to true then other models and parameters will not update correctly.

        See Also
        --------
        OrbitBase.set_orbit
        """

        planet_pointer = self.find_planet_pointer(planet_reference)
        planet_loc = planet_pointer.orbit_location

        if planet_pointer is self.host and not self.star_host and not set_by_planet:
            if debug_mode:
                # This would change the orbit between the host and the star, not between the target body and the star.
                # This may not be what the user wants to do.
                log("Attempting to change the host planet's orbit", level='debug')

        new_semi_major_axis = value_cleanup(new_semi_major_axis)

        if semi_major_axis_in_au:
            new_semi_major_axis = Au2m(new_semi_major_axis)

        self._semi_major_axis[planet_loc] = new_semi_major_axis
        # Changing the orbital semi-major axis also changes the orbital frequency. Update via Kepler Laws
        if planet_pointer is self.host:
            self._orbital_freqs[planet_loc] = \
                np.sqrt(G * (planet_pointer.mass + self.star.mass) / new_semi_major_axis**3)
        else:
            # Star and Target bodies all orbit the host (albeit the star doesn't do much...).
            self._orbital_freqs[planet_loc] = \
                np.sqrt(G * (planet_pointer.mass + self.host.mass) / new_semi_major_axis**3)
        if not set_by_planet:
            self.update_orbit(planet_reference, set_by_planet=False)

    def set_eccentricity(self, planet_reference: PlanetRefType, new_eccentricity: FloatArray,
                         set_by_planet: bool = False):
        """ Set the eccentricity of a planet at planet_loc

        Use Orbit.set_orbit if setting more than one state parameter.

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference used to find the planet
        new_eccentricity : FloatArray
            Eccentricity (only elliptical and circular orbits are supported in TidalPy)
            Optional
        set_by_planet : bool = False
            A flag used by TidalPy planets to make calls to this method. Leave it False for user use.
            If set to true then other models and parameters will not update correctly.

        See Also
        --------
        OrbitBase.set_orbit
        """

        planet_pointer = self.find_planet_pointer(planet_reference)
        planet_loc = planet_pointer.orbit_location

        if planet_pointer is self.host and not self.star_host and not set_by_planet:
            if debug_mode:
                # This would change the orbit between the host and the star, not between the target body and the star.
                # This may not be what the user wants to do.
                log("Attempting to change the host planet's orbit", level='debug')

        new_eccentricity = value_cleanup(new_eccentricity)

        self._eccentricities[planet_loc] = new_eccentricity
        if not set_by_planet:
            self.update_orbit(planet_reference, set_by_planet=False)

    def set_inclination(self, planet_reference: PlanetRefType, new_inclination: FloatArray,
                        inclination_in_deg: bool = False, set_by_planet: bool = False):
        """ Set the inclination of a planet at planet_loc

        Use Orbit.set_orbit if setting more than one state parameter.

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference used to find the planet
        new_inclination : FloatArray
            Orbital inclination relative to the orbital plane in [rads]
            Optional
        inclination_in_deg : bool = False
            If true then TidalPy will convert the input into radians (the unit that is required for calculations)
        set_by_planet : bool = False
            A flag used by TidalPy planets to make calls to this method. Leave it False for user use.
            If set to true then other models and parameters will not update correctly.

        See Also
        --------
        OrbitBase.set_orbit
        """

        planet_pointer = self.find_planet_pointer(planet_reference)
        planet_loc = planet_pointer.orbit_location

        if planet_pointer is self.host and not self.star_host and not set_by_planet:
            if debug_mode:
                # This would change the orbit between the host and the star, not between the target body and the star.
                # This may not be what the user wants to do.
                log("Attempting to change the host planet's orbit", level='debug')

        new_inclination = value_cleanup(new_inclination)

        if inclination_in_deg:
            new_inclination = np.deg2rad(new_inclination)

        self._inclinations[planet_loc] = new_inclination
        if not planet_reference:
            self.update_orbit(planet_reference, set_by_planet=False)

    # Functionality to get the orbit state variables
    def get_orbital_freq(self, planet_reference: PlanetRefType) -> np.ndarray:

        planet_loc = self.find_planet_pointer(planet_reference).orbit_location

        if planet_loc == 1:
            planet_loc = self.host_tide_raiser_loc

        return self._orbital_freqs[planet_loc]

    def get_semi_major_axis(self, planet_reference: PlanetRefType) -> np.ndarray:

        planet_loc = self.find_planet_pointer(planet_reference).orbit_location

        if planet_loc == 1:
            planet_loc = self.host_tide_raiser_loc

        return self._semi_major_axis[planet_loc]

    def get_eccentricity(self, planet_reference: PlanetRefType) -> np.ndarray:

        planet_loc = self.find_planet_pointer(planet_reference).orbit_location

        if planet_loc == 1:
            planet_loc = self.host_tide_raiser_loc

        return self._eccentricities[planet_loc]

    def get_inclination(self, planet_reference: PlanetRefType) -> np.ndarray:

        planet_loc = self.find_planet_pointer(planet_reference).orbit_location

        if planet_loc == 1:
            planet_loc = self.host_tide_raiser_loc

        return self._inclinations[planet_loc]

    def get_derivative_eccentricity(self, planet_reference: PlanetRefType) -> np.ndarray:

        if not self.is_time_study:
            raise IncorrectModelInitialized('To access time derivatives the orbit must be initialized as a time domain '
                                            'study.')
        planet_loc = self.find_planet_pointer(planet_reference).orbit_location
        if planet_loc == 0:
            raise ParameterValueError()
        elif planet_loc == 1:
            planet_loc = self.host_tide_raiser_loc

        return self._derivative_e[planet_loc]

    def get_derivative_semi_major(self, planet_reference: PlanetRefType) -> np.ndarray:

        if not self.is_time_study:
            raise IncorrectModelInitialized('To access time derivatives the orbit must be initialized as a time domain '
                                            'study.')
        planet_loc = self.find_planet_pointer(planet_reference).orbit_location
        if planet_loc == 0:
            raise ParameterValueError()
        elif planet_loc == 1:
            planet_loc = self.host_tide_raiser_loc

        return self._derivative_a[planet_loc]

    def get_derivative_spin(self, planet_reference: PlanetRefType) -> np.ndarray:

        planet_ref = self.find_planet_pointer(planet_reference)
        return planet_ref.derivative_spin

    def update_orbit(self, planet_reference: PlanetRefType = None, set_by_planet: bool = False):
        """ Updates the state variables of the orbit instance

        Whenever changes to the orbit or to its planet's are made then the orbit's state variables may need updating.
        This method will perform all necessary updates.

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference to the planet which has changed. If set to None then update_orbit will be applied to all
            tidal objects
        set_by_planet : bool
            Was this method called by a planet instance
        """

        if planet_reference is None:
            planet_pointers = self.target_bodies
        else:
            planet_pointers = [self.find_planet_pointer(planet_reference)]

        last_p_index = len(planet_pointers) - 1
        for p_index, planet_pointer in enumerate(planet_pointers):
            planet_loc = planet_pointer.orbit_location
            if not set_by_planet:
                # Need to tell the planet to update any orbital-dependent state
                planet_pointer.update_orbit()

                if self.has_tidal_host and p_index == last_p_index:
                    if planet_loc == self.host_tide_raiser_loc:
                        # The host needs to be updated too (only once all the other planets have been updated).
                        self.host.update_orbit()

            # Update time derivatives
            if self.is_time_study:
                target_tidal_heating = planet_pointer.tidal_heating
                target_tidal_torque = planet_pointer.tidal_ztorque
                target_spin_freq = planet_pointer.spin_freq
                semi_major_axis = self._semi_major_axis[planet_loc]
                eccentricity = self._eccentricities[planet_loc]
                target_mass = planet_pointer.mass
                host_mass = self.host.mass

                if target_tidal_heating is None:
                    raise ParameterMissingError
                if target_tidal_torque is None:
                    raise ParameterMissingError
                if target_spin_freq is None:
                    raise ParameterMissingError

                if self.duel_dissipation:
                    host_tidal_heating = self.host.tidal_heating
                    host_tidal_torque = self.host.tidal_ztorque
                    host_spin_freq = self.host.spin_freq

                    if host_tidal_heating is None:
                        raise ParameterMissingError
                    if host_tidal_torque is None:
                        raise ParameterMissingError
                    if host_spin_freq is None:
                        raise ParameterMissingError

                    da_dt_func = diff_eqs_duel_dissipation['semi_major_axis']
                    de_dt_func = diff_eqs_duel_dissipation['eccentricity']
                    self._derivative_a[planet_loc] = \
                        da_dt_func(semi_major_axis, host_mass, target_mass,
                                   target_spin_freq, target_tidal_torque, target_tidal_heating,
                                   host_spin_freq, host_tidal_torque, host_tidal_heating)
                    self._derivative_e[planet_loc] = \
                        de_dt_func(semi_major_axis, eccentricity, host_mass, target_mass,
                                   target_spin_freq, target_tidal_torque, target_tidal_heating,
                                   host_spin_freq, host_tidal_torque, host_tidal_heating)
                else:
                    da_dt_func = diff_eqs_single_dissipation['semi_major_axis']
                    de_dt_func = diff_eqs_single_dissipation['eccentricity']
                    self._derivative_a[planet_loc] = \
                        da_dt_func(semi_major_axis, host_mass, target_mass,
                                   target_spin_freq, target_tidal_torque, target_tidal_heating)
                    self._derivative_e[planet_loc] = \
                        de_dt_func(semi_major_axis, eccentricity, host_mass, target_mass,
                                   target_spin_freq, target_tidal_torque, target_tidal_heating)

    def calculate_insolation(self, planet_reference: PlanetRefType, set_planet_param: bool = True):
        """ Calculate the insolation heating received by a planet located at planet_loc

        The star-planet separation (semi-major axis) and eccentricity are used to estimate the orbit-averaged
        insolation heating received at the surface of the target planet. If the star is the host of the system then
        the target body's semi-a and eccentricity will be used. Otherwise the orbit.host's parameters will be used.

        The actual method used to make the calculation is stored in the target body's equilibrium_insolation_func.
        It is set by the planet's configuration---the difference between the models is how they handle an eccentric
        orbit.

        Parameters
        ----------
        planet_reference : PlanetRefType
            Reference used to find the planet
        set_planet_param : bool = True
            If the method should set the planet's insolation_heating state variable or not

        Returns
        -------
        insolation_heating : FloatArray
            The orbit averaged heating received at the surface of the target planet in [Watts]
        """

        planet_pointer = self.find_planet_pointer(planet_reference)
        planet_loc = planet_pointer.orbit_location

        if planet_pointer is self.star:
            raise ParameterValueError('Can not calculation insolation heating for the star.')

        # These separations and eccentricities should be relative to the star!
        star_separation = self.equilib_distance_funcs[planet_loc]()
        if star_separation is None:
            raise ParameterMissingError

        star_eccentricity = self.equilib_eccentricity_funcs[planet_loc]()
        if star_eccentricity is None:
            if debug_mode:
                log('Attempting to calculate insolation heating with no eccentricity set. Using e=0.')
            star_eccentricity = np.asarray(0.)

        insolation_heating = planet_pointer.equilibrium_insolation_func(
                self.star.luminosity, star_separation, planet_pointer.albedo, planet_pointer.radius,
                star_eccentricity)

        if set_planet_param:
            planet_pointer.insolation_heating = insolation_heating

        return insolation_heating

    @property
    def host_tide_raiser_loc(self) -> int:
        return self._host_tide_raiser_loc

    @host_tide_raiser_loc.setter
    def host_tide_raiser_loc(self, value: int):

        if value > len(self.all_objects):
            raise ParameterValueError('Host tide raiser location must be the orbit location of one of the target bodies.')
        elif value in [0, 1]:
            raise ParameterValueError('Host tide raiser location can not be the host or stars locations')
        elif value < 0:
            raise ParameterValueError('Host tide raiser location must be positive')

        self._host_tide_raiser_loc = value

        # Now update the host planet with relevant information
        pseudo_host_reference = self.all_objects[self.host_tide_raiser_loc]
        self.host.tidal_susceptibility_inflated = (3. / 2.) * G * pseudo_host_reference.mass**2 * \
                                                  self.host.radius**5
        self.host.tide_raiser_ref = pseudo_host_reference

    def __iter__(self):
        return iter(self.all_objects)

    @staticmethod
    def rads2days(radians_per_second: FloatArray) -> FloatArray:
        """ Convert radians/sec (frequency) to days (period)

        Wrapper for TidalPy.utilities.conversions.rads2days
        """

        return rads2days(radians_per_second)

    @staticmethod
    def days2rads(days: FloatArray) -> FloatArray:
        """ Convert days (period) to radians/sec (frequency)

        Wrapper for TidalPy.utilities.conversions.days2rads
        """

        return days2rads(days)
