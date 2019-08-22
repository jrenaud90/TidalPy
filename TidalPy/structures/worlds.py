from __future__ import annotations

import copy
import os
from typing import List, TYPE_CHECKING, Tuple

import burnman
import dill
import numpy as np

from ..dynamics.modes import nsr_modes, spin_sync_modes
from ..utilities.numpy_help import value_cleanup
from .defaults import world_defaults
from .physical import PhysicalObjSpherical
from .. import debug_mode, use_disk, tidalpy_dir, __version__
from ..configurations import (auto_save_planet_config_to_rundir, auto_save_planet_config_to_tidalpydir,
                              auto_save_planet_dill_to_rundir, auto_save_planet_dill_to_tidalpydir, overwrite_configs,
                              overwrite_dills)
from ..dynamics import spin_rate_derivative
from ..exceptions import (ImproperAttributeHandling, ParameterMissingError, ReinitError, UnusualRealValueError)
from ..graphics.planet_plot import geotherm_plot
from ..initialize import log
from ..io import inner_save_dir
from ..physics.stellar import (efftemp_from_luminosity, equilibrium_insolation_functions, equilibrium_temperature,
                               luminosity_from_efftemp, luminosity_from_mass)
from ..utilities.dill_helper import dill_file_path
from ..structures.layers import ThermalLayer

planet_config_loc = os.path.join(tidalpy_dir, 'planets', 'planet_configs')

if TYPE_CHECKING:
    from ..orbit.orbit import OrbitBase


class WorldBase(PhysicalObjSpherical):
    class_type = 'base'
    default_config = copy.deepcopy(world_defaults)

    def __init__(self, world_config: dict, call_reinit: bool = True):

        # Load in defaults
        self.default_config = self.default_config[self.class_type]

        super().__init__(config=world_config, call_reinit=False)

        # Pull out switch information
        self.is_spin_sync = self.config['force_spin_sync']

        # Independent State variables
        self._spin_freq = None
        self._time = None

        # Dependent State Variables
        self._insolation_heating = None
        self._surface_temperature = None

        # Orbit reference
        self._orbit = None  # type: OrbitBase or None

        # Other flags
        self.is_host = False

        # Additional orbit information
        # orbit_location is an index for where a planet is at in the orbit's geometry. This is set by the initialization
        # of an orbit class.
        #   0 == Star location
        #   1 == Host location
        #   2+ == Various target bodies
        self.orbit_location = None

        # Models to be initialized later
        self.equilibrium_insolation_func = None

        # Constants to be set in reinit
        self.emissivity = None
        self.albedo = None
        self._name = None
        self.pyname = f'{self.name}_{self.class_type}'

        if call_reinit:
            self.reinit()

    def reinit(self):

        super().reinit()

        self.name = self.config['name']
        log(f'Reinit called for planet: {self.name}', level='debug')

        # Pull out constants
        self.emissivity = self.config['emissivity']
        self.albedo = self.config['albedo']

        # Setup functions and models
        self.equilibrium_insolation_func = equilibrium_insolation_functions[self.config['equilibrium_insolation_model']]

    def set_geometry(self, radius: float, mass: float):

        super().set_geometry(radius, mass, thickness=radius)

    def update_orbit(self):
        """ Performs state updates whenever the planet's orbital parameters are changed
        """

        # Tell the orbit class to update this planet's insolation heating
        self.orbit.calculate_insolation(self.orbit_location)

    def set_state(self, orbital_freq: np.ndarray = None, semi_major_axis: np.ndarray = None,
                  eccentricity: np.ndarray = None, inclination: np.ndarray = None,
                  spin_freq: np.ndarray = None, time: np.ndarray = None):
        """ Set multiple orbital parameters at once, this reduces the number of calls to self.orbit_change

        This contains a wrapper to the orbit class method set_state. It extends that function by including the spin
        frequency.

        This function has better performance than the individual setters for these parameters if (and only if) you are
        changing two or more of them at the same time.

        Parameters
        ----------
        orbital_freq : np.ndarray
            Mean orbital frequency of an object in [rads s-1]
            Optional, exclusive w/ semi_major_axis
        semi_major_axis : np.ndarray
            Semi-major axis of an object in [m]
            Optional, exclusive w/ semi_major_axis
        eccentricity : np.ndarray
            Orbital eccentricity
            Optional
        inclination : np.ndarray
            Orbital inclination in [rads]
            Optional
        spin_freq : np.ndarray
            Object's rotation frequency in [rads s-1]
            Optional
        time : np.ndarray
            Object's calculation time (used in things like radiogenics) in [s]
            Optional
        """

        if time is not None:
            self._time = value_cleanup(time)

        self.orbit.set_orbit(self.orbit_location, orbital_freq, semi_major_axis, eccentricity, inclination,
                             set_by_planet=True)
        if spin_freq is not None:
            if type(spin_freq) != np.ndarray:
                spin_freq = np.asarray([spin_freq])
            self._spin_freq = spin_freq

        self.update_orbit()

    def kill_world(self):
        """ Performs saving tasks when the world is about to be deleted due to end of run

        The exit_planets variable in the main TidalPy configurations controls rather or not this method ever gets
        called automatically.
        """

        log(f'Killing world {self.name}...', level='debug')
        # Save configuration file
        if use_disk:
            if auto_save_planet_config_to_tidalpydir:
                tidalpy_planet_cfg_dir = [planet_config_loc]
            else:
                tidalpy_planet_cfg_dir = list()
            self.save_config(save_to_run_dir=auto_save_planet_config_to_rundir,
                             additional_save_dirs=tidalpy_planet_cfg_dir, overwrite=overwrite_configs)

        # Dill this object and save it
        if use_disk:
            dill_name, tidalpy_dill_path = dill_file_path(self.name)
            if auto_save_planet_dill_to_tidalpydir:
                if os.path.isfile(tidalpy_dill_path) and not overwrite_dills:
                    pass
                else:
                    with open(tidalpy_dill_path, 'wb') as dill_file:
                        dill.dump(self, dill_file)
            if auto_save_planet_dill_to_rundir:
                rundir_dill_path = os.path.join(inner_save_dir, dill_name)
                if os.path.isfile(rundir_dill_path) and not overwrite_dills:
                    pass
                else:
                    with open(rundir_dill_path, 'wb') as dill_file:
                        dill.dump(self, dill_file)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value: str):
        assert type(value) == str
        self._name = value
        self.pyname = f'{self.name}_{self.class_type}'

    @property
    def orbit(self) -> OrbitBase:
        return self._orbit

    @orbit.setter
    def orbit(self, value):
        raise ImproperAttributeHandling('Orbit is automatically passed to a planet when a new orbit class is created.')

    # Independent state variables
    @property
    def time(self) -> np.ndarray:
        return self._time

    @time.setter
    def time(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])

        if debug_mode:
            if np.any(abs(value)) > 1.e6:
                raise UnusualRealValueError(f'Time should be entered in units of [Myr]. |{value}| seems very large.')

        self._time = value_cleanup(value)

    @property
    def spin_freq(self) -> np.ndarray:
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])
        self._spin_freq = value_cleanup(value)
        if not self.is_spin_sync and self.orbit is not None:
            # NSR will change tidal modes
            self.update_orbit()

    @property
    def insolation_heating(self) -> np.ndarray:
        return self._insolation_heating

    @insolation_heating.setter
    def insolation_heating(self, value: np.ndarray):
        self._insolation_heating = value

    @property
    def surface_temperature(self) -> np.ndarray:
        return self._surface_temperature

    @surface_temperature.setter
    def surface_temperature(self, value: np.ndarray):
        raise ImproperAttributeHandling

    # Wrappers for the orbit class's state variable setters and getters
    @property
    def semi_major_axis(self) -> np.ndarray:
        return self.orbit.get_semi_major_axis(self.orbit_location)

    @semi_major_axis.setter
    def semi_major_axis(self, value: np.ndarray):
        self.orbit.set_semi_major_axis(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def orbital_freq(self) -> np.ndarray:
        return self.orbit.get_orbital_freq(self.orbit_location)

    @orbital_freq.setter
    def orbital_freq(self, value: np.ndarray):

        self.orbit.set_orbital_freq(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def eccentricity(self) -> np.ndarray:
        return self.orbit.get_eccentricity(self.orbit_location)

    @eccentricity.setter
    def eccentricity(self, value: np.ndarray):
        self.orbit.set_eccentricity(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    @property
    def inclination(self) -> np.ndarray:
        return self.orbit.get_inclination(self.orbit_location)

    @inclination.setter
    def inclination(self, value: np.ndarray):
        self.orbit.set_inclination(self.orbit_location, value, set_by_planet=True)
        self.update_orbit()

    # Aliased Parameter Names
    @property
    def orbital_motion(self):
        return self.orbital_freq

    @orbital_motion.setter
    def orbital_motion(self, value):
        self.orbital_freq = value

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'


class BasicWorld(WorldBase):
    class_type = 'basic'

    def __init__(self, planet_config: dict, call_reinit: bool = True):
        super().__init__(planet_config, call_reinit=call_reinit)

        self.set_geometry(planet_config['radius'], planet_config['mass'])

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'


class GasGiant(BasicWorld):
    class_type = 'gasgiant'

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'


class Star(BasicWorld):
    class_type = 'star'

    def __init__(self, star_config: dict, call_reinit: bool = True):

        super().__init__(star_config, call_reinit=False)

        # Star-specific attributes
        self.luminosity = None
        self.effective_temperature = None

        if call_reinit:
            self.reinit()

    def reinit(self):

        super().reinit()

        self.luminosity = self.config.get('luminosity', None)
        self.effective_temperature = self.config.get('effective_temperature', None)

        if self.luminosity is None:
            # If no luminosity provided: Try to convert effective surface temperature
            if self.effective_temperature is None:
                # if that fails, try to estimate from mass
                log(
                        'Luminosity and effective temperature of {self.name} was not provided. Estimating from stellar mass.')
                self.luminosity = luminosity_from_mass(self.mass)
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)
            else:
                self.luminosity = luminosity_from_efftemp(self.effective_temperature, self.radius)
        else:
            if self.effective_temperature is None:
                self.effective_temperature = efftemp_from_luminosity(self.luminosity, self.radius)

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'


class TidalWorld(WorldBase):
    class_type = 'tidal'

    def __init__(self, planet_config: dict, burnman_world: burnman.Planet, bm_layers: list, call_reinit: bool = True):

        super().__init__(planet_config, call_reinit=False)

        # Store burnman information
        self.bm_world = burnman_world
        self.bm_layers = bm_layers

        # Setup geometry based on burnman
        self.set_geometry(self.bm_world.radius_planet, self.bm_world.mass)

        # Pull out other burnman information
        self._moi = self.bm_world.moment_of_inertia

        # Parameters only applicable to a tidally-active host
        self.tide_raiser_ref = None
        self.use_real_moi = None
        self.fixed_q = None

        # Setup layers
        self.layers_byname = dict()
        self.layers = list()
        for layer_i, bm_layer in enumerate(self.bm_layers):
            layer_name = bm_layer.name
            layer_config = self.config['layers'][layer_name]
            layer = ThermalLayer(layer_name, self, bm_layer, layer_config)
            if layer_i == len(self.bm_layers) - 1:
                # Top most layer
                layer.is_top_layer = True
            self.layers.append(layer)
            # Store both the title and lower case name of the layer (these have the same pointer)
            self.layers_byname[layer_name.title()] = layer
            self.layers_byname[layer_name.lower()] = layer
            setattr(self, layer.name.title(), layer)
            setattr(self, layer.name.lower(), layer)

        # Constants used in calculations
        self.tidal_susceptibility_inflated = None  # = 1.5 * G * M_pri^2 * R_sec^5 -- Set by orbit class

        # Provide references to other layers
        for layer_below, layer, layer_above in zip([None] + self.layers[:-1], self.layers, self.layers[1:] + [None]):
            layer.layer_below = layer_below
            layer.layer_above = layer_above

        # Concatenate the slice information from the layers
        gravity_list = [layer.gravity_slices for layer in self]
        pressure_list = [layer.pressure_slices for layer in self]
        density_list = [layer.density_slices for layer in self]
        radius_list = [layer.radii for layer in self]
        self.gravity_slices = np.concatenate(gravity_list)
        self.pressure_slices = np.concatenate(pressure_list)
        self.density_slices = np.concatenate(density_list)
        self.radii = np.concatenate(radius_list)

        # Dependent state variables
        self._tidal_modes = None
        self._tidal_freqs = None
        self._tidal_heating_coeffs = None
        self._tidal_ztorque_coeffs = None
        self._global_love = None
        self._tidal_heating = None
        self._tidal_ztorque = None
        self._tidal_susceptibility = None
        self._derivative_spin = None

        if call_reinit:
            self.reinit()

    def reinit(self, update_spin: bool = True, update_orbit: bool = True):

        super().reinit()

        # Save the TidalPy version for used to make the world.
        self._config['TidalPy_version'] = __version__

        # Load in some global parameters
        if update_spin:
            spin_freq = self.config.get('spin_freq', None)
            spin_period = self.config.get('spin_period', None)
            if spin_freq is not None and spin_period is not None:
                log(f'Both spin frequency and period were provided for {self.name}. '
                    f'Using frequency instead.', level='info')
            if spin_freq is None and spin_period is not None:
                # Assume spin period is in days
                spin_freq = 2. * np.pi / (spin_period * 24. * 60. * 60.)
            self._spin_freq = np.asarray([spin_freq])

        # Load in model flags
        self.use_real_moi = self.config['use_real_moi']

        # Load in other parameters
        self.fixed_q = self.config['quality_factor']

        # Check configs
        if self._old_config is None:
            reload = False
        else:
            reload = True

        if reload:
            # Determine if planet must be made again due to config changes that impact burnman results
            for layer_name, layer_dict in self.config['layers'].items():
                old_layer_dict = self._old_config['layers'][layer_name]
                for critical_attribute in ['material', 'material_source', 'type']:
                    if layer_dict[critical_attribute] != old_layer_dict[critical_attribute]:
                        raise ReinitError

        # Load in new configs into layers and call their reinits
        for layer in self:
            if layer.config is not None:
                new_config = {**layer.config, **self.config['layers'][layer.name]}
            else:
                new_config = self.config['layers'][layer.name]
            layer.user_config = new_config
            layer.update_config()
            layer.reinit()

            # Store the layer config in the world so that you can see the defaults that are loaded in
            self.config['layers'][layer.name] = layer.config

        # Try to initialize as much as we can
        for layer in self:
            layer.update_tides()

        if self.orbit is not None and update_orbit:
            self.orbit.update_orbit(self, set_by_planet=True)

    def find_layer(self, layer_name: str) -> ThermalLayer:
        """ Returns a reference to a layer with the provided name

        Layers are also stored in the planet's __dict__ and can be accessed via planet.<layer_name>

        Parameters
        ----------
        layer_name : str
            Name assigned to layer in the planet's original configuration

        Returns
        -------
        layer : ThermalLayer
            Reference to the layer class
        """

        layer = self.layers_byname[layer_name]
        return layer

    def find_layer_byradius(self, radius: float) -> ThermalLayer:
        """ Returns a reference to a layer that the provided radius resides in

        If the provided radius is at the interface of two layers this method will choose the lower layer.

        Parameters
        ----------
        radius : float
            Radius in [m] at which to search for a layer

        Returns
        -------
        layer : ThermalLayer
            Reference to the layer class at that radius

        """

        assert radius <= self.radius

        for layer in self.layers:
            if layer.radius >= radius:
                return layer
        raise LookupError()

    def update_orbit(self):
        """ Updates the planet's state when there has been a change to its orbit

        Whenever the planet's orbital frequency, eccentricity, inclination, and/or spin frequency change this method
        should be called to update other models and properties.

        Currently it updates the tidal modes and their corresponding coefficients used in the calculation of
        tidal heating and tidal torques.

        See Also
        --------
        TidalWorld.update_global_tides
        """

        if self.is_host:
            # Orbital parameters for a host planet are pulled from its tide raising target body (set by the orbit class)
            orbit_reference_object = self.tide_raiser_ref
        else:
            orbit_reference_object = self

        # Pull out orbit information
        eccentricity = orbit_reference_object.eccentricity
        inclination = orbit_reference_object.inclination
        orbital_freq = orbit_reference_object.orbital_freq
        semi_major_axis = orbit_reference_object.semi_major_axis

        # Spin rate is stored in self regardless if the planet is a host or not
        spin_freq = self.spin_freq

        # Preform checks
        # These are not required
        if eccentricity is None:
            eccentricity = np.asarray(0.)
        if inclination is None:
            inclination = np.asarray(0.)
        # These are required
        if orbital_freq is None:
            raise ParameterMissingError
        if semi_major_axis is None:
            raise ParameterMissingError

        if self.is_spin_sync:
            self._tidal_modes, self._tidal_freqs, self._tidal_heating_coeffs, self._tidal_ztorque_coeffs = \
                spin_sync_modes(orbital_freq, eccentricity, inclination)
        else:
            if self.spin_freq is None:
                raise ParameterMissingError

            self._tidal_modes, self._tidal_freqs, self._tidal_heating_coeffs, self._tidal_ztorque_coeffs = \
                nsr_modes(orbital_freq, spin_freq, eccentricity, inclination)

        # The semi-major axis should have changed. We can now update the tidal susceptibility
        self._tidal_susceptibility = self.tidal_susceptibility_inflated / semi_major_axis**6

        # Now update calculations at the layer level
        for layer in self:
            layer.update_tides(call_from_world=True)

        # Perform an update to the global tides as well since layer love numbers have likely changed
        self.update_global_tides()

        super().update_orbit()

    def update_global_tides(self, set_layer_heating: bool = True):
        """ Updates the planet's tidal state whenever there is a change to the orbit or to a layer's thermal state

        Changes to a planet's orbit, or to its constituent layers' thermal state, will change the global love number and
        tidal heating & torque. It is easier for a planet to make these calculations and pass on the relevant results to
        its layers.

        The layer's complex love number is stored as a Tuple[np.ndarray] where their is an item stored in the tuple for
        each tidal mode. For a planet that is forced to be spin-sync'd, there will only be one member. NSR worlds can
        have many tidal modes. These tidal modes will be collapsed down to one number within this method. The collapse
        procedure requires multiplying the Love number by a coefficient which is unique for each mode. This process is
        done with the list comprehensions half-way through this method.

        The global Love number is simply a sum of the mode-collapsed, layer Love numbers.

        Anytime a layer Love number is used to find a global value it is first multiple by a scaling factor called the
        'Tidal Volume Fraction' (denoted by the variable 'scale' below). Scale is used to mimic the affect of a
        multilayer calculation in a 1D model. By default it is a volume scale of the layer's tidally active region
        compared to the total planet's volume. This can be turned off in the configuration by setting
        'use_tvf' to False.

        Parameters
        ----------
        set_layer_heating : bool
            If true then this method will automatically propagate tidal heating to the relevant layers.

        See Also
        --------
        TidalWorld.update_orbit
        """

        # The layer love numbers should have already been updated before this point, so we can simply pull their values
        layer_loves = [layer.complex_love for layer in self]  # type: List[Tuple[np.ndarray]]

        global_love_list = list()
        global_heating_list = list()
        global_ztorque_list = list()

        for layer, love_by_mode in zip(self, layer_loves):
            scale = layer.tidal_scale
            if layer.is_tidal:
                layer_heating_bymodes = [-np.imag(love) * heating_coeff for love, heating_coeff in
                                         zip(love_by_mode, self.tidal_heating_coeffs)]
                layer_ztorque_bymodes = [-np.imag(love) * ztorque_coeff for love, ztorque_coeff in
                                         zip(love_by_mode, self.tidal_ztorque_coeffs)]
                layer_heating_sum = self.tidal_susceptibility * scale * sum(layer_heating_bymodes)
                layer_ztorque_sum = self.tidal_susceptibility * scale * sum(layer_ztorque_bymodes)
            else:
                layer_ztorque_sum = 0.
                layer_heating_sum = 0.

            # Global Values
            global_love_list.append(scale * sum(love_by_mode))
            global_heating_list.append(layer_heating_sum)
            global_ztorque_list.append(layer_ztorque_sum)

            # Update layer values
            if set_layer_heating:
                layer.tidal_heating = layer_heating_sum

        self._global_love = sum(global_love_list)
        self._tidal_heating = sum(global_heating_list)
        self._tidal_ztorque = sum(global_ztorque_list)

        # Update Spin Derivative
        if self.use_real_moi:
            moi = self.moi
            if moi is None:
                log('Tried to use real MOI but it was not set: using ideal instead.', level='debug')
                moi = self.moi_ideal
        else:
            moi = self.moi_ideal
        if moi is None:
            raise ParameterMissingError
        self._derivative_spin = spin_rate_derivative(self.tidal_ztorque, moi)

    def update_surface_temperature(self):
        """

        Returns
        -------
        surface_temperature : FloatArray
            Equilibrium surface temperature in [K]
        """

        # A change in the insolation heating will change the surface temperature which may impact cooling
        heat_flux_out = self.layers[-1].heat_flux
        if heat_flux_out is None:
            heat_flux_out = np.asarray([0.])
        heating_out = heat_flux_out * self.surface_area_outer
        total_heating = heating_out + self.insolation_heating
        surface_temperature = equilibrium_temperature(total_heating, self.radius, self.emissivity)

        # TODO: this may need to be a convergence since this propagation will change the heat flux and thus the
        #    equilibrium temp.
        #    Right now I am leaving the below commented out because update_cooling makes a call to this setter in the
        #    top layer. If this is uncommented then there will be circular calls
        #    This is probably safe to leave out for now since the surface temp is going to be changing so little
        #    due to changes in the internal temp.
        # for layer in self[::-1]:
        #     # Propagate the changes downward through the planet
        #     layer.update_cooling()

        self._surface_temperature = surface_temperature

        return surface_temperature

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

    @property
    def time(self) -> np.ndarray:
        return self._time

    @time.setter
    def time(self, value: np.ndarray):
        if type(value) != np.ndarray:
            value = np.asarray([value])

        if debug_mode:
            if np.any(abs(value)) > 1.e6:
                raise UnusualRealValueError(f'Time should be entered in units of [Myr]. |{value}| seems very large.')

        self._time = value

        # Update radiogenics
        for layer in self:
            if layer.radiogenics is not None:
                layer.update_radiogenics(force_calculation=True)

    @property
    def insolation_heating(self) -> np.ndarray:
        return self._insolation_heating

    @insolation_heating.setter
    def insolation_heating(self, value: np.ndarray):
        self._insolation_heating = value
        self.update_surface_temperature()

    @property
    def tidal_modes(self):
        return self._tidal_modes

    @tidal_modes.setter
    def tidal_modes(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_freqs(self):
        return self._tidal_freqs

    @tidal_freqs.setter
    def tidal_freqs(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating_coeffs(self):
        return self._tidal_heating_coeffs

    @tidal_heating_coeffs.setter
    def tidal_heating_coeffs(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_ztorque_coeffs(self):
        return self._tidal_ztorque_coeffs

    @tidal_ztorque_coeffs.setter
    def tidal_ztorque_coeffs(self, value):
        raise ImproperAttributeHandling

    @property
    def global_love(self):
        return self._global_love

    @global_love.setter
    def global_love(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_heating(self):
        return self._tidal_heating

    @tidal_heating.setter
    def tidal_heating(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_ztorque(self):
        return self._tidal_ztorque

    @tidal_ztorque.setter
    def tidal_ztorque(self, value):
        raise ImproperAttributeHandling

    @property
    def tidal_susceptibility(self):
        return self._tidal_susceptibility

    @tidal_susceptibility.setter
    def tidal_susceptibility(self, value):
        raise ImproperAttributeHandling

    @property
    def derivative_spin(self):
        return self._derivative_spin

    @derivative_spin.setter
    def derivative_spin(self, value):
        raise ImproperAttributeHandling

    def __iter__(self):
        """ Planet will iterate over its layers

        Returns
        -------
        iter(self.layers)
            The iterator of the layer list.
        """

        return iter(self.layers)

    def __repr__(self):

        if 'name' in self.__dict__:
            if self.name is not None:
                return f'{self.name} {self.__class__} object at {hex(id(self))}'
        return f'{self.__class__} object at {hex(id(self))}'


world_types = {
    'star'    : Star,
    'basic'   : BasicWorld,
    'simple'  : BasicWorld,
    'gas'     : GasGiant,
    'gasgiant': GasGiant,
    'burnman' : TidalWorld,
    'tidal'   : TidalWorld
}
