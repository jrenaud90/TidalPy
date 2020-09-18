from __future__ import annotations

import os
from typing import TYPE_CHECKING, Union

import numpy as np

from .defaults import world_defaults
from .. import PhysicalObjSpherical
from ... import debug_mode, use_disk, tidalpy_loc, configurations, log
from ...exceptions import (ImproperPropertyHandling, UnusualRealValueError,
                           AttributeNotSetError, IOException, ConfigPropertyChangeError,
                           IncorrectMethodToSetStateProperty, UnknownModelError, InitiatedPropertyChangeError)
from ...helpers.orbit_help import pull_out_orbit_from_config
from ...stellar import equilibrium_insolation_functions, EquilibFuncType, calc_equilibrium_temperature
from ...utilities.graphics import geotherm_plot
from ...utilities.types import FloatArray, NoneType

planet_config_loc = os.path.join(tidalpy_loc, 'planets', 'planet_configs')

if TYPE_CHECKING:
    from ...orbit import Orbit


class BaseWorld(PhysicalObjSpherical):

    """ WorldBase Class - Base class used to build other world classes.


    See Also
    --------
    Parent Class:
        TidalPy.structures.PhysicalObjSpherical
    Child Classes:
        TidalPy.structures.worlds.TidalWorld
        TidalPy.structures.worlds.GasGiantWorld
        TidalPy.structures.worlds.StarWorld
        TidalPy.structures.worlds.LayeredWorld
        TidalPy.structures.worlds.GasGiantLayeredWorld
        TidalPy.structures.worlds.BurnManWorld
    """

    default_config = world_defaults
    world_class = 'base'

    def __init__(self, world_config: dict, name: str = None, initialize: bool = True):
        """ BaseWorld constructor

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

        # Load in defaults
        self.default_config = self.default_config[self.world_class]

        # Key Attributes
        self._name = name
        log.debug(f'Setting up new world: {self.name}; class type = {self.world_class}.')

        super().__init__(config=world_config)

        # Configuration variables
        self._albedo = None
        self._emissivity = None
        self._force_spin_sync = False
        self._equilibrium_insolation_func = None
        self._internal_to_surf_heating_frac = None

        # Independent state variables
        self._spin_frequency = None
        self._obliquity = None
        self._time = None
        self._surf_temperature = None
        self._internal_to_surf_heating = None

        # Orbit reference
        self.orbit = None  # type: Union[NoneType, Orbit]
        self.tidal_host = None

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False, reinit_geometry: bool = True, set_by_burnman: bool = False):
        """ Initialize or Reinitialize the world based on changes to its configurations.

        This must be called at least once before an instance can be used. The constructor will automatically make an
            initial call to reinit unless told to not to.

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this method has been called (additional steps may be
                preformed during the first reinit call).
        reinit_geometry : bool = True
            If `True`, the initializer will automatically call the `set_geometry()` method.
        set_by_burnman : bool = False
            Set to `True` if called from a burnman world.
        """

        super().reinit(initial_init=initial_init, set_by_burnman=set_by_burnman)

        # Pull out some basic information from the config
        if self.name is None and 'name' in self.config:
            self.name = self.config['name']

        # Set flags
        self._force_spin_sync = self.config['force_spin_sync']

        # Set Thermal configurations
        self._albedo = self.config['albedo']
        self._emissivity = self.config['emissivity']

        # Set Physical configurations
        self._pressure_above = self.config.get('surface_pressure', 0.)

        # Set Orbit-related configurations
        self._obliquity = self.config.get('obliquity', 0.)

        # Set surface temperature / stellar interaction configurations
        self._internal_to_surf_heating_frac = self.config['fraction_internal_heating_to_surface']
        insol_equilib_func = self.config['equilibrium_insolation_model']
        try:
            self._equilibrium_insolation_func = equilibrium_insolation_functions[insol_equilib_func]
        except KeyError:
            log.error(f'Unknown equilibrium insolation function model encountered in {self}.')
            raise UnknownModelError

        # Setup geometry
        if reinit_geometry:
            self.set_geometry(self.config['radius'], self.config['mass'])
            self.set_static_pressure(pressure_above=self.pressure_above, build_slices=True)

        # Clean up config:
        if 'radii' in self.config:
            del self._config['radii']

        if self.orbit is not None:
            # Update orbit with any new configurations
            orbital_freq, semi_major_axis, eccentricity = pull_out_orbit_from_config(self.config)
            self.orbit.set_state(self, new_eccentricity=eccentricity,
                                 new_orbital_frequency=orbital_freq, new_semi_major_axis=semi_major_axis,
                                 called_from_orbit=True)


            # Update orbit with config changes.
            self.update_orbit()

    def set_geometry(self, radius: float, mass: float, thickness: float = None, mass_below: float = 0.,
                     update_state_geometry: bool = True, build_slices: bool = True):
        """ Calculates and sets the world's physical parameters based on user provided input.

        Assumptions
        -----------
        Spherical Geometry

        Parameters
        ----------
        radius : float
            Outer radius of object [m]
        mass : float
            Mass of object [kg]
        thickness : float = None
            Thickness of the object [m]
        mass_below : float = 0.
            Mass below this object (only applicable for shell-like structures)
            Used in gravity and pressure calculations
        update_state_geometry : bool = True
            Update the class' state geometry
        build_slices : bool = True
            If True, method will attempt to calculate gravities, densities, etc. for each slice.

        """

        # Thickness of a world will always be equal to its radius.
        del thickness
        super().set_geometry(radius, mass, thickness=radius, mass_below=0.,
                             update_state_geometry=update_state_geometry, build_slices=build_slices)

    def update_time(self):
        """ Update any attributes or methods that rely on the current time.

        This can be called by a Orbit class if one is set.
        """

        log.debug(f'Update time method called for {self}.')

        # Most update_time functionality implemented by child classes.

    def update_surface_temperature(self):
        """ Update the planet's surface temperature based on its distance from its host star.

        Orbit must be applied before surface temperature can be updated.
        """

        log.debug(f'Update surface temperature method called for {self}.')

        if self.orbit is not None:
            self._surf_temperature = \
                calc_equilibrium_temperature(self.insolation_heating, self.radius, self.internal_to_surf_heating,
                                             self.emissivity, self.internal_to_surf_heating_frac)

    def update_orbit(self):
        """ Performs state updates whenever the planet's orbital parameters are changed
        """

        log.debug(f'Update orbit method called for {self}.')

        # Tell the orbit class to update this planet's insolation heating
        if self.orbit is not None:
            self.orbit.calculate_insolation(self.orbit_location)

        # A change to semi-major axis for a planet orbiting a star will cause a change to its surface temperature
        self.update_surface_temperature()

    def update_tides(self):
        """ Update the tides module

        Only applicable for some child class types
        """

        log.debug(f'Update tides method called for {self}.')

        # Most update_tides functionality implemented by child classes.

    def clear_state(self, preserve_orbit: bool = False):
        """ Clear the world's current state variables back to their defaults.

        The defaults may be Nones or set by the user-provided configuration.

        Parameters
        ----------
        preserve_orbit: bool = False
            If `True`, data about this planet's orbit will be cleared from any associated Orbit classes.
        """

        # The parent class only does a log, since we are doing that here there is no need to call parent method.
        log.debug(f'Clear state called for {self}. Orbit preserved = {preserve_orbit}.')

        # Clear world-specific state properties
        self._spin_frequency = None
        self._obliquity = None
        self._time = None
        self._surf_temperature = None
        self._internal_to_surf_heating = None


        # Setup functions and models
        self.equilibrium_insolation_func = equilibrium_insolation_functions[self.config['equilibrium_insolation_model']]

        # Reset any orbits this world is connected to
        if not preserve_orbit and self.orbit is not None:
            # Purge orbital properties for this world from the Orbit class.
            self.orbit.clear_state(clear_all=False, clear_specific=self, clear_world_state=False)

    def set_state(self, spin_frequency: FloatArray = None, obliquity: FloatArray = None, time: FloatArray = None,
                  orbital_frequency: FloatArray = None, orbital_period: FloatArray = None,
                  semi_major_axis: FloatArray = None, eccentricity: FloatArray = None, set_by_world: bool = False):
        """ Set multiple orbital parameters at once, this reduces the number of calls to self.orbit_change

        This contains a wrapper to the orbit class method set_state. It extends that function by including the spin
        frequency.

        This function has better performance than the individual setters for these parameters if (and only if) you are
        changing two or more of them at the same time.

        Parameters
        ----------
        spin_frequency : FloatArray = None
            New spin frequency for the world [rad s-1].
        obliquity : FloatArray = None
            New obliquity for the world relative to its orbit around the tidal host [rad].
        time : FloatArray = None
            Time used in integration and radiogenic calculations [Myr].
        orbital_frequency : FloatArray = None
            New orbital frequency (orbital motion) for the world around its tidal host [rad s-1].
        orbital_period : FloatArray = None
            New orbital period for the world around its tidal host [days].
        semi_major_axis: FloatArray = None
            New orbital separation between world and its tidal host [m].
        eccentricity: FloatArray = None
            New orbital eccentricity relative to the tidal host.
        set_by_world: bool = False
            If `True`, then worlds update methods will not be called.
        """

        # Track if we need to update tides or not.
        update_tides_flag = False
        orbit_will_call_update_tides_flag = False
        update_time_flag = False

        # Make flags
        new_spin_frequency = spin_frequency is not None
        new_obliquity = obliquity is not None
        new_time = time is not None
        new_orbital_frequency = orbital_frequency is not None
        new_orbital_period = orbital_period is not None
        new_semi_major_axis = semi_major_axis is not None
        new_eccentricity = eccentricity is not None

        # Check for time changes
        if new_time:
            self.set_time(time, call_updates=False)
            update_time_flag = True

        # Check for world properties that can effect tides but not the orbit directly
        if new_spin_frequency:
            self.set_spin_frequency(spin_frequency, call_updates=False)
            update_tides_flag = True

        if new_obliquity:
            self.set_obliquity(obliquity, call_updates=True)
            update_tides_flag = True

        # Check if any orbital updates
        if any((new_orbital_frequency, new_orbital_period, new_semi_major_axis, new_eccentricity)):
            # Tell orbit to update the state of this world.
            self.orbit.set_state(self, new_eccentricity=eccentricity, new_semi_major_axis=semi_major_axis,
                                 new_orbital_frequency=orbital_frequency, new_orbital_period=orbital_period)
            update_tides_flag = True
            orbit_will_call_update_tides_flag = True

        if not set_by_world:
            # Check what updates need to be called, if any.
            if update_time_flag:
                self.update_time()

            if update_tides_flag:
                # Tides need to be updated, but is the orbit going to call the method anyways?
                if not orbit_will_call_update_tides_flag:
                    self.update_tides()

    def set_time(self, new_time: FloatArray, call_updates: bool = True):
        """ Set the time of the world.

        Parameters
        ----------
        new_time : FloatArray
            Time used in integration and radiogenic calculations [Myr]
        call_updates : bool = True
            If `True`, method will call the update time method.
        """

        if self.orbit is None:
            if debug_mode:
                if np.any(abs(new_time) > 1.e6):
                    raise UnusualRealValueError(f'Time should be entered in units of [Myr]. '
                                                f'|{new_time}| seems very large.')

            self._time = new_time
            if not call_updates:
                self.update_time()
        else:
            raise ImproperPropertyHandling('Time must be set at the Orbit-level once an orbit is applied.')

    def set_spin_frequency(self, new_spin_frequency: FloatArray, call_updates: bool = True):
        """ Update the world's spin frequency.

        Parameters
        ----------
        new_spin_frequency : FloatArray
            New spin frequency for the world [rad s-1]
        call_updates : bool = True
            If `True`, method will call the update tides method.
        """

        if debug_mode:
            if np.any(np.abs(new_spin_frequency) > 1.e-3):
                raise UnusualRealValueError(f'Spin-frequency should be entered in units of [rads s-1]. '
                                            f'|{new_spin_frequency}| seems very large.')

        self._spin_frequency = new_spin_frequency

        if call_updates:
            self.update_tides()

    def set_obliquity(self, new_obliquity: FloatArray, call_updates: bool = True):
        """ Set the world's obliquity.

        This obliquity must be relative to the orbital plane defined by the tidal target and tidal host [1]_.
            If the star is neither the host nor target, then it should not be used as a reference object.

        References
        ----------
        .. [1] J. P. Renaud et al, "Tidal Dissipation in Dual-Body, Highly Eccentric, and
           Non-synchronously Rotating Systems: Applications to Pluto-Charon and the Exoplanet TRAPPIST-1e"
           The Planetary Science Journal, vol. 22, pp. TBA, 2020.

        Parameters
        ----------
        new_obliquity : FloatArray
            New obliquity for the world relative to its orbit around the tidal host [rad]
        call_updates : bool = True
            If `True`, method will not call the update tides method.
        """

        if debug_mode:
            if np.any(new_obliquity > 7.):
                raise UnusualRealValueError('Obliquity should be entered in radians. '
                                            f'A value of {np.max(new_obliquity)} seems unusually large.')

        self._obliquity = new_obliquity

        if not call_updates:
            self.update_tides()

    def paint(self, depth_plot: bool = False, auto_show: bool = True):
        """ Create a geotherm or depth plot of the planet's gravity, pressure, and density
        Parameters
        ----------
        depth_plot : bool = False
            If true the plot will be versus depth rather than radius.
        auto_show : bool = False
            Calls plt.show() if true.
        Returns
        -------
        figure: matplotlib.pyplot.figure
        """

        figure = geotherm_plot(self.radii, self.gravity_slices, self.pressure_slices, self.density_slices,
                               bulk_density=self.density_bulk, planet_name=self.name,
                               planet_radius=self.radius, depth_plot=depth_plot, auto_show=auto_show)

        return figure

    def save_world(self, save_dir: str = None, no_cwd: bool = False, save_to_tidalpy_dir: bool = False):
        """ Save the world's configuration file to a specified directory.

        Parameters
        ----------
        save_dir : str = None
            The current working directory will be preappended unless no_cwd is set to true.
            If no directory is provided it will be saved os.getcwd()
        no_cwd : bool = False
            If True, the current working directory will *not* be prepended to the provided directory.
        save_to_tidalpy_dir : bool = False
            If True, the config will be saved to the TidalPy directory as well as the CWD.
        """

        if not use_disk:
            raise IOException("User attempted to save a world's configurations when TidalPy has been set to "
                              "not use_disk")

        log.debug(f'Saving world configurations for {self}.')

        save_locales = list()
        if save_to_tidalpy_dir:
            save_locales.append(planet_config_loc)

        if save_dir is not None:
            if no_cwd:
                save_locales.append(save_dir)
            else:
                save_locales.append(os.path.join(os.getcwd(), save_dir))
        else:
            save_locales.append(os.getcwd())

        # No need to save to run dir as that will automatically happen when the planet is killed.
        self.save_config(save_to_run_dir=False, additional_save_dirs=save_locales,
                         overwrite=configurations['overwrite_configs'])

    def kill_world(self):
        """ Performs saving tasks when the world is about to be deleted due to end of run

        The exit_planets variable in the main TidalPy configurations controls rather or not this method ever gets
        called automatically.
        """

        log.debug(f'Killing world {self}.')

        # Save configuration file
        if use_disk:
            if configurations['auto_save_planet_config_to_tidalpydir']:
                tidalpy_planet_cfg_dir = [planet_config_loc]
            else:
                tidalpy_planet_cfg_dir = list()
            self.save_config(save_to_run_dir=configurations['auto_save_planet_config_to_rundir'],
                             additional_save_dirs=tidalpy_planet_cfg_dir,
                             overwrite=configurations['overwrite_configs'])



    # # Initialized properties
    @property
    def name(self) -> str:
        """ Name of the world """
        return self._name

    @name.setter
    def name(self, value):
        raise InitiatedPropertyChangeError


    # # Configuration properties
    @property
    def albedo(self) -> float:
        """ World's Albedo """
        return self._albedo

    @albedo.setter
    def albedo(self, value):
        raise ConfigPropertyChangeError

    @property
    def emissivity(self) -> float:
        """ World's Emissivity """
        return self._emissivity

    @emissivity.setter
    def emissivity(self, value):
        raise ConfigPropertyChangeError

    @property
    def force_spin_sync(self) -> bool:
        """ Flag that is used to force the world's spin rate to equal its orbital motion if changed. """
        return self._force_spin_sync

    @force_spin_sync.setter
    def force_spin_sync(self, value):
        raise ConfigPropertyChangeError

    @property
    def equilibrium_insolation_func(self) -> EquilibFuncType:
        """ Flag that is used to force the world's spin rate to equal its orbital motion if changed. """
        return self._equilibrium_insolation_func

    @equilibrium_insolation_func.setter
    def equilibrium_insolation_func(self, value):
        raise ConfigPropertyChangeError

    @property
    def internal_to_surf_heating_frac(self) -> EquilibFuncType:
        """ Fraction of internal heating that makes its way to the surface (used for surface equilibrium temperature
        calculations). """
        return self._internal_to_surf_heating_frac

    @internal_to_surf_heating_frac.setter
    def internal_to_surf_heating_frac(self, value):
        raise ConfigPropertyChangeError


    # # State properties
    @property
    def time(self) -> FloatArray:
        """ The time of either the Orbit or the object (used for radiogenic calculations) [Myr] """
        if self.orbit is None:
            return self._time
        else:
            return self.orbit.universal_time

    @time.setter
    def time(self, new_time):
        self.set_time(new_time)

    @property
    def spin_frequency(self) -> FloatArray:
        """ Spin (Sidereal Rotation) Frequency of the World [rad s-1] """
        return self._spin_frequency

    @spin_frequency.setter
    def spin_frequency(self, new_spin_frequency: FloatArray):
        self.set_spin_frequency(new_spin_frequency)

    @property
    def surface_temperature(self) -> FloatArray:
        """ World's Surface Temperature [K] """
        # TODO: Alias this with PhysicalObject's self.temperature_outer?
        return self._surf_temperature

    @surface_temperature.setter
    def surface_temperature(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def obliquity(self) -> FloatArray:
        """ World's Obliquity [rad]

        This obliquity must be relative to the orbital plane defined by the tidal target and tidal host [1]_.
            If the star is neither the host nor target, then it should not be used as a reference object.

        References
        ----------
        .. [1] J. P. Renaud et al, "Tidal Dissipation in Dual-Body, Highly Eccentric, and
           Non-synchronously Rotating Systems: Applications to Pluto-Charon and the Exoplanet TRAPPIST-1e"
           The Planetary Science Journal, vol. 22, pp. TBA, 2020.
        """

        return self._obliquity

    @obliquity.setter
    def obliquity(self, new_obliquity: FloatArray):
        self.set_obliquity(new_obliquity)

    @property
    def internal_to_surf_heating(self) -> np.ndarray:
        return self._internal_to_surf_heating

    @internal_to_surf_heating.setter
    def internal_to_surf_heating(self, value):
        raise IncorrectMethodToSetStateProperty


    # Outer-scope properties
    # # Orbit Class
    @property
    def semi_major_axis(self):
        """ World's Semi-Major Axis (stored in the world's Orbit class) [m] """
        return self.orbit.get_semi_major_axis(self)

    @semi_major_axis.setter
    def semi_major_axis(self, new_semi_major_axis: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set semi-major axis until an Orbit class has been applied to {self}.')

        self.orbit.set_semi_major_axis(self, new_semi_major_axis)

    @property
    def orbital_frequency(self):
        """ World's Orbital Frequency (inverse of orbital period; stored in the world's Orbit class) [rad s-1]"""
        return self.orbit.get_orbital_frequency(self)

    @orbital_frequency.setter
    def orbital_frequency(self, new_orbital_frequency: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set orbital frequency until an Orbit class has been applied to {self}.')

        self.orbit.set_orbital_frequency(self, new_orbital_frequency)

    @property
    def orbital_period(self):
        """ World's Orbital Period (inverse of orbital frequency; stored in the world's Orbit class) [days]"""
        return self.orbit.get_orbital_period(self)

    @orbital_period.setter
    def orbital_period(self, new_orbital_period: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set orbital period until an Orbit class has been applied to {self}.')

        self.orbit.set_orbital_period(self, new_orbital_period)

    @property
    def eccentricity(self):
        """ World's Orbital Eccentricity (inverse of orbital period; stored in the world's Orbit class) """
        return self.orbit.get_eccentricity(self)

    @eccentricity.setter
    def eccentricity(self, new_eccentricity: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set orbital eccentricity until an Orbit class has been applied '
                                       f'to {self}.')

        self.orbit.set_eccentricity(self, new_eccentricity)

    @property
    def insolation_heating(self):
        """ World's Surface Heating due to Stellar Insolation [W] """
        if self.orbit is not None:
            return self.orbit.get_insolation(self)
        return None

    @insolation_heating.setter
    def insolation_heating(self, value):
        raise IncorrectMethodToSetStateProperty


    # # Aliased properties
    @property
    def spin_freq(self):
        """ Alias of BaseWorld.spin_frequency """
        return self.spin_frequency

    @spin_freq.setter
    def spin_freq(self, value):
        self.spin_frequency = value

    @property
    def orbital_freq(self):
        """ Alias of BaseWorld.orbital_frequency """
        return self.orbital_frequency

    @orbital_freq.setter
    def orbital_freq(self, new_orbital_frequency):
        self.orbital_frequency = new_orbital_frequency

    @property
    def n(self):
        """ Alias of BaseWorld.orbital_frequency """
        return self.orbital_frequency

    @n.setter
    def n(self, new_orbital_frequency):
        self.orbital_frequency = new_orbital_frequency

    @property
    def orbital_motion(self):
        """ Alias of BaseWorld.orbital_frequency """
        return self.orbital_frequency

    @orbital_motion.setter
    def orbital_motion(self, value):
        self.orbital_frequency = value


    # Dunder properties
    def __str__(self):

        name = self.name
        if name is None:
            name = f'Unknown World'

        name = f'{name} ({self.world_class} world)'

        return name

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name}, {self.__class__} object at {hex(id(self))}'
        return f'Unknown World, {self.__class__} object at {hex(id(self))}'
