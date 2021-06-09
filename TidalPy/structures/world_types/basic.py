from __future__ import annotations

import os
from typing import TYPE_CHECKING, Union

import numpy as np

from .defaults import world_defaults
from .. import PhysicalObjSpherical
from ..helpers.orbit_config import pull_out_orbit_from_config
from ... import debug_mode, use_disk, tidalpy_loc, configurations, log
from ...exceptions import (ImproperPropertyHandling, UnusualRealValueError,
                           AttributeNotSetError, IOException, ConfigPropertyChangeError,
                           IncorrectMethodToSetStateProperty, UnknownModelError, InitiatedPropertyChangeError,
                           OuterscopePropertySetError)
from ...stellar import calc_equilibrium_temperature, equilibrium_insolation_functions, EquilibFuncType
from ...toolbox.conversions import days2rads, rads2days
from ...utilities.graphics import geotherm_plot
from ...utilities.types import FloatArray, NoneType

planet_config_loc = os.path.join(tidalpy_loc, 'planets', 'planet_configs')

if TYPE_CHECKING:
    from . import AllWorldType
    from ..orbit import Orbit


class BaseWorld(PhysicalObjSpherical):

    """ WorldBase Class - Base class used to build other world methods.


    See Also
    --------
    Parent Class:
        TidalPy.structures.PhysicalObjSpherical
    Child Classes:
        TidalPy.structures.world_types.TidalWorld
        TidalPy.structures.world_types.GasGiantWorld
        TidalPy.structures.world_types.StarWorld
        TidalPy.structures.world_types.LayeredWorld
        TidalPy.structures.world_types.GasGiantLayeredWorld
        TidalPy.structures.world_types.BurnManWorld
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
            Please see files stored in <TidalPy directory>/structures/world_configs for example configuration dict.
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
        self._spin_period = None
        self._obliquity = None
        self._time = None
        self._surface_temperature = None
        self._insolation_heating = None

        # Orbit reference
        self.orbit = None  # type: Union[NoneType, Orbit]

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
            # TODO: Isn't this done in the orbit class?
            orbital_freq, semi_major_axis, eccentricity = pull_out_orbit_from_config(self.config)
            self.orbit.set_state(self, eccentricity=eccentricity,
                                 orbital_frequency=orbital_freq, semi_major_axis=semi_major_axis,
                                 call_orbit_change=False)

    def update_surface_temperature(self, called_from_cooling: bool = False):
        """ Surface temperature has changed - Perform any calculations that may have also changed.

        Parameters
        ----------
        called_from_cooling : bool = False
            Flag to avoid recursive loops between surface temperature and cooling.
        """

        log.debug(f'Updated surface temperature method called for {self}.')

        internal_heating = self.get_internal_heating_to_surface()

        if self.insolation_heating is not None:
            self._surface_temperature = \
                calc_equilibrium_temperature(self.insolation_heating, self.radius,
                                             internal_heating, self.emissivity)

        self.surface_temperature_changed(called_from_cooling=called_from_cooling)

    def time_changed(self):
        """ The world's time has been changed. Make any necessary updates.
        """

        log.debug(f'Time changed for {self}.')

        # Most updated_time functionality implemented by child methods.

    def orbit_spin_changed(self, orbital_freq_changed: bool = False, spin_freq_changed: bool = False,
                           eccentricity_changed: bool = False, obliquity_changed: bool = False,
                           call_orbit_dissipation: bool = True):
        """ The world's orbit, spin, and/or obliquity has changed. Make any necessary updates """

        log.debug(f'Orbit, spin, and/or obliquity changed called for {self}.')

        # More updates implemented by child classes.

    def surface_temperature_changed(self, called_from_cooling: bool = False):
        """ Surface temperature has changed - Perform any calculations that may have also changed.

        Parameters
        ----------
        called_from_cooling : bool = False
            Flag to avoid recursive loops between surface temperature and cooling.
        """

        log.debug(f'Surface temperature changed for {self}.')

        # More updates implemented by child classes.

    def clear_state(self, preserve_orbit: bool = False):
        """ Clear the world's current state variables back to their defaults.

        The defaults may be Nones or set by the user-provided configuration.

        Parameters
        ----------
        preserve_orbit: bool = False
            If `True`, data about this planet's orbit will be cleared from any associated Orbit methods.
        """

        # The parent class only does a log, since we are doing that here there is no need to call parent method.
        log.debug(f'Clear state called for {self}. Orbit preserved = {preserve_orbit}.')

        # Clear world-specific state properties
        self._spin_frequency = None
        self._obliquity = None
        self._time = None
        self._surface_temperature = None

        # Setup functions and models
        self._equilibrium_insolation_func = equilibrium_insolation_functions[self.config['equilibrium_insolation_model']]

        # Reset any orbits this world is connected to
        if not preserve_orbit and self.orbit is not None:
            # Purge orbital properties for this world from the Orbit class.
            self.orbit.clear_state(clear_all=False, clear_specific=self, clear_world_state=False)

    def set_state(self, spin_frequency: FloatArray = None, spin_period: FloatArray = None,
                  obliquity: FloatArray = None, time: FloatArray = None,
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
        spin_period : FloatArray = None
            New spin period for the world [days]
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
            If `True`, then world_types update methods will not be called.
        """

        log.debug(f'Set state called for {self}.')

        # Make flags
        new_spin_frequency = spin_frequency is not None
        new_spin_period = spin_period is not None
        new_obliquity = obliquity is not None
        new_time = time is not None
        new_orbital_frequency = orbital_frequency is not None
        new_orbital_period = orbital_period is not None
        new_semi_major_axis = semi_major_axis is not None
        new_eccentricity = eccentricity is not None

        # Check for time changes
        if new_time:
            self.set_time(time, call_updates=False)

        # Check for world properties that can effect tides but not the orbit directly
        if new_spin_frequency:
            if new_spin_period:
                log.warning(f'Both a new spin frequency and a spin period or provided to {self}. Using spin frequency.')
            self.set_spin_frequency(spin_frequency, call_updates=False)
        elif new_spin_period:
            self.set_spin_period(spin_period, call_updates=False)
            # A change to the spin period will also change the spin frequency
            new_spin_frequency = True

        if new_obliquity:
            self.set_obliquity(obliquity, call_updates=False)

        # Check if any orbital updates
        if any((new_orbital_frequency, new_orbital_period, new_semi_major_axis, new_eccentricity)):
            # Tell orbit to update the state of this world.
            self.orbit.set_state(self, eccentricity=eccentricity, semi_major_axis=semi_major_axis,
                                 orbital_frequency=orbital_frequency, orbital_period=orbital_period,
                                 set_by_world=True)

            if new_orbital_period and not new_orbital_frequency:
                # A change to the orbital period will also change the spin frequency
                new_orbital_frequency = True
            if new_semi_major_axis and not new_orbital_frequency:
                # A change to the orbital semi-major axis will also change the spin frequency
                new_orbital_frequency = True

        if not set_by_world:

            if new_orbital_frequency and self.force_spin_sync:
                # If the orbital frequency changes and the world is forced into synchronous rotation then its spin
                #    frequency needs to be updated.
                self.set_spin_frequency(self.orbital_frequency, call_updates=False)

            if new_time:
                self.time_changed()

            if any((new_spin_frequency, new_spin_period, new_orbital_frequency, new_obliquity, new_eccentricity)):
                self.orbit_spin_changed(spin_freq_changed=new_spin_frequency,
                                        orbital_freq_changed=new_orbital_frequency,
                                        obliquity_changed=new_obliquity,
                                        eccentricity_changed=new_eccentricity)

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

    def set_time(self, time: FloatArray, call_updates: bool = True):
        """ Set the time of the world.

        Parameters
        ----------
        time : FloatArray
            Time used in integration and radiogenic calculations [Myr]
        call_updates : bool = True
            If `True`, method will call the update time method.
        """

        if self.orbit is None:
            self._time = time
            if call_updates:
                self.time_changed()
        else:
            raise ImproperPropertyHandling('Time must be set at the Orbit-level once an orbit is applied.')

    def set_spin_frequency(self, spin_frequency: FloatArray, call_updates: bool = True):
        """ Update the world's spin frequency.

        Parameters
        ----------
        spin_frequency : FloatArray
            New spin frequency for the world [rad s-1]
        call_updates : bool = True
            If `True`, method will call the update tides method.
        """

        if debug_mode:
            if np.any(np.abs(spin_frequency) > 1.e-3):
                raise UnusualRealValueError(f'Spin-frequency should be entered in units of [rads s-1]. '
                                            f'|{spin_frequency}| seems very large.')

        self._spin_frequency = spin_frequency
        self._spin_period = rads2days(spin_frequency)
        if call_updates:
            self.orbit_spin_changed(spin_freq_changed=True)

    def set_spin_period(self, spin_period: FloatArray, call_updates: bool = True):
        """ Update the world's spin period in days.

        Parameters
        ----------
        spin_period : FloatArray
            New spin period for the world [days]
        call_updates : bool = True
            If `True`, method will call the update tides method.
        """

        spin_frequency = days2rads(spin_period)

        self.set_spin_frequency(spin_frequency, call_updates=call_updates)

    def set_obliquity(self, obliquity: FloatArray, call_updates: bool = True):
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
        obliquity : FloatArray
            New obliquity for the world relative to its orbit around the tidal host [rad]
        call_updates : bool = True
            If `True`, method will not call the update tides method.
        """

        if debug_mode:
            if np.any(obliquity > 7.):
                raise UnusualRealValueError('Obliquity should be entered in radians. '
                                            f'A value of {np.max(obliquity)} seems unusually large.')

        self._obliquity = obliquity
        if call_updates:
            self.orbit_spin_changed(obliquity_changed=True)

    def set_insolation_heating(self, insolation_heating: FloatArray):
        """ Set a new insolation heating received on the world's surface from a host star.

        This method will update the surface temperature which will in turn change the thermal state.

        See Also
        --------
        TidalPy.orbit.physics.PhysicsOrbit.calculate_insolation()

        Parameters
        ----------
        insolation_heating : FloatArray
            Heating received from a stellar host [W]
        """

        self._insolation_heating = insolation_heating

        # Now recalculate the surface equilibrium temperature
        self.update_surface_temperature()

    def get_internal_heating_to_surface(self) -> Union[NoneType, FloatArray]:
        """ Get the amount of internal heating that is making it to the surface.

        Returns
        -------
        internal_heating_to_surface : Union[NoneType, FloatArray]
            Amount of heating that is reaching the surface [W]
        """

        internal_heating_to_surface = None

        # Child methods override this with more functionality

        return internal_heating_to_surface

    def paint(self, depth_plot: bool = False, auto_show: bool = True, return_fig: bool = False):
        """ Create a geotherm or depth plot of the planet's gravity, pressure, and density
        Parameters
        ----------
        depth_plot : bool = False
            If `True` the plot will be versus depth rather than radius.
        auto_show : bool = False
            Calls plt.show() if true.
        return_fig : bool = False
            If `True`, return the matplotlib fig object otherwise return True.

        Returns
        -------
        figure: matplotlib.pyplot.figure
        """

        figure = geotherm_plot(self.radii, self.gravity_slices, self.pressure_slices, self.density_slices,
                               bulk_density=self.density_bulk, planet_name=self.name,
                               planet_radius=self.radius, depth_plot=depth_plot, auto_show=auto_show)

        if return_fig:
            return figure
        else:
            return True

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
    def spin_period(self) -> FloatArray:
        """ Spin (Sidereal Rotation) Period of the World [days] """
        return self._spin_period

    @spin_period.setter
    def spin_period(self, new_spin_period: FloatArray):
        self.set_spin_period(new_spin_period)

    @property
    def surface_temperature(self) -> FloatArray:
        """ World's Surface Temperature [K] """
        return self._surface_temperature

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
    def insolation_heating(self) -> FloatArray:
        """ Surface heating received on a world from its host star [W] """
        return self._insolation_heating

    @insolation_heating.setter
    def insolation_heating(self, new_insolation_heating: FloatArray):
        self.set_insolation_heating(new_insolation_heating)


    # Outer-scope properties
    # # Orbit Class
    @property
    def tidal_host(self) -> Union[NoneType, 'AllWorldType']:
        """ Wrapper for the orbit class's get_tidal_host method """
        if self.orbit is None:
            return None
        else:
            return self.orbit.get_tidal_host(self)

    @tidal_host.setter
    def tidal_host(self, value):
        raise OuterscopePropertySetError

    @property
    def semi_major_axis(self):
        """ World's Semi-Major Axis (stored in the world's Orbit class) [m] """
        if self.orbit is None:
            return None
        return self.orbit.get_semi_major_axis(self)

    @semi_major_axis.setter
    def semi_major_axis(self, new_semi_major_axis: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set semi-major axis until an Orbit class has been applied to {self}.')

        self.orbit.set_semi_major_axis(self, new_semi_major_axis)

    @property
    def orbital_frequency(self):
        """ World's Orbital Frequency (inverse of orbital period; stored in the world's Orbit class) [rad s-1]"""
        if self.orbit is None:
            return None
        return self.orbit.get_orbital_frequency(self)

    @orbital_frequency.setter
    def orbital_frequency(self, new_orbital_frequency: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set orbital frequency until an Orbit class has been applied to {self}.')

        self.orbit.set_orbital_frequency(self, new_orbital_frequency)

    @property
    def orbital_period(self):
        """ World's Orbital Period (inverse of orbital frequency; stored in the world's Orbit class) [days]"""
        if self.orbit is None:
            return None
        return self.orbit.get_orbital_period(self)

    @orbital_period.setter
    def orbital_period(self, new_orbital_period: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set orbital period until an Orbit class has been applied to {self}.')

        self.orbit.set_orbital_period(self, new_orbital_period)

    @property
    def eccentricity(self):
        """ World's Orbital Eccentricity (stored in the world's Orbit class) """
        if self.orbit is None:
            return None
        return self.orbit.get_eccentricity(self)

    @eccentricity.setter
    def eccentricity(self, new_eccentricity: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set orbital eccentricity until an Orbit class has been applied '
                                       f'to {self}.')

        self.orbit.set_eccentricity(self, new_eccentricity)

    @property
    def stellar_eccentricity(self):
        """ World's Orbital Eccentricity relative to its host star (stored in the world's Orbit class) """
        if self.orbit is None:
            return None
        return self.orbit.get_stellar_eccentricity(self)

    @stellar_eccentricity.setter
    def stellar_eccentricity(self, new_eccentricity: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set the stellar eccentricity until an Orbit class has been applied '
                                       f'to {self}.')
        self.orbit.set_stellar_eccentricity(self, new_eccentricity)

    @property
    def stellar_distance(self):
        """ World's Orbital Distance relative to its host star (stored in the world's Orbit class) [m] """
        if self.orbit is None:
            return None
        return self.orbit.get_stellar_distance(self)

    @stellar_distance.setter
    def stellar_distance(self, new_distance: FloatArray):

        if self.orbit is None:
            raise AttributeNotSetError(f'Can not set the stellar distance until an Orbit class has been applied '
                                       f'to {self}.')
        self.orbit.set_stellar_distance(self, new_distance)

    @property
    def eccentricity_time_derivative(self) -> Union[NoneType, FloatArray]:
        """ Derivative of eccentricity with respect to time (only effects due to tides are considered) """
        if self.orbit is None:
            return None
        else:
            return self.orbit.get_eccentricity_time_derivative(self)

    @eccentricity_time_derivative.setter
    def eccentricity_time_derivative(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def semi_major_axis_time_derivative(self) -> Union[NoneType, FloatArray]:
        """ Derivative of the semi-major axis with respect to time (only effects due to tides are considered) """
        if self.orbit is None:
            return None
        else:
            return self.orbit.get_semi_major_axis_time_derivative(self)

    @semi_major_axis_time_derivative.setter
    def semi_major_axis_time_derivative(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def orbital_motion_time_derivative(self) -> Union[NoneType, FloatArray]:
        """ Derivative of the orbital mean motion with respect to time (only effects due to tides are considered) """
        if self.orbit is None:
            return None
        else:
            return self.orbit.get_orbital_motion_time_derivative(self)


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

        str_ = f'{name} ({self.__class__.__name__})'

        return str_