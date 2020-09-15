from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

from .basic import LayerBase
from ...cooling import CoolingModel
from ...exceptions import (ImproperPropertyHandling, ParameterMissingError, OuterscopePropertySetError)
from ...radiogenics.radiogenics import Radiogenics
from ...rheology import Rheology
from ...utilities.numpyHelper.array_shape import reshape_help
from ...utilities.types import NoneType, FloatArray

if TYPE_CHECKING:
    from ..worlds import TidalWorldType


class PhysicsLayer(LayerBase):

    """ PhysicsLayer
    Layer object to store parameters geometric and physical properties calculated by TidalPy based on a user-provided
        configuration dictionary. PhysicsLayer contains additional properties and methods used to perform various
        thermal and tidal calculations.

    Notes:
    .. Should not be used with Burnman calculated layers (see BurnmanLayer instead)

    See Also
    --------
    TidalPy.structures.layers.BurnmanLayer
    """

    layer_class = 'physics'

    def __init__(self, layer_name: str, layer_index: int, world: 'TidalWorldType', layer_config: dict,
                 is_top_layer: bool, initialize: bool = True):
        """ Physics layer constructor

        Parameters
        ----------
        layer_name : str
            User-friendly name of layer.
        layer_index : int
            Location of layer within a world (0 indicates center-most).
        world : LayeredWorldType
            World instance where layer was initialized in.
        layer_config : dict
            Layer's user-provided configurations.
        is_top_layer : bool
            If `True`, this layer is the top-most layer.
        initialize : bool = True
            If `True`, then the Layer's reinit is called at the end of the constructor.
        """

        # Setup Physical Layer Geometry
        super().__init__(layer_name, layer_index, world, layer_config, is_top_layer, initialize=False)

        # State Derivatives
        self._deriv_temperature = None

        # Other attributes
        self.use_pressure_in_strength_calc = True

        # Model holders - setup in reinit()
        self.rheology = None       # type: Union['Rheology', NoneType]
        self.cooling_model = None  # type: Union['CoolingModel', NoneType]
        self.radiogenics = None    # type: Union['Radiogenics', NoneType]

        # Set up a pressure that will persist if the layer's state is cleared
        # TODO: should this be reinit?
        self._persistent_pressure = self.pressure

        # Material Properties
        #     These are calculated by BurnMan in the BurnmanLayer version, can be set by user/config here.
        self.bulk_modulus = None
        self.thermal_expansion = None
        self.specific_heat = None
        self.energy_per_therm = None

        # Material properties that might have been affected by new configuration file (BurnMan does not calculate these)
        self.static_shear_modulus = None
        self.thermal_diffusivity = None
        self.thermal_conductivity = None
        self.heat_fusion = None
        self.temp_ratio = None
        self.stefan = None

        if initialize:
            self.reinit(initial_init=initialize)

    def reinit(self, initial_init: bool = False, set_by_burnman: bool = False, initialize_geometry: bool = True):
        """ Reinitialize the physical object by pulling in any potentially new configurations

        Parameters
        ----------
        initial_init : bool = False
            Set to `True` for the first time an instance is created.
        set_by_burnman : bool = False
            Set to `True` if a Burnman layer/world constructor is calling reinit
        initialize_geometry : bool = False
            Set to `True` if the set_geometry method should be called from within reinit
        """

        # Base class's reinit is called *after* the geometry is set (so that tidal volume fraction is set correctly)
        super().reinit(initial_init=initial_init, set_by_burnman=set_by_burnman,
                       initialize_geometry=initialize_geometry)

        # Material properties that might have been affected by new configuration files
        self.static_shear_modulus = self.config['shear_modulus']
        self.heat_fusion = self.config['heat_fusion']
        self.thermal_conductivity = self.config['thermal_conductivity']
        self.temp_ratio = self.config['boundary_temperature_ratio']
        self.stefan = self.config['stefan']

        # Setup Models
        self.rheology = Rheology(self, store_config_in_layer=True)
        self.radiogenics = Radiogenics(self, store_config_in_layer=True)
        self.cooling_model = CoolingModel(self, store_config_in_layer=True)

        # Find all heating sources
        self.heat_sources = list()
        if 'off' not in self.radiogenics.model:
            self.heat_sources.append(lambda : self.radiogenic_heating)
        if self.is_tidal:
            self.heat_sources.append(lambda: self.tidal_heating)

    def clear_state(self, clear_pressure: bool = False):

        super().clear_state(clear_pressure=clear_pressure)

        # State Derivatives
        self._deriv_temperature = None

        # Clear the state of all inner-scope classes
        for model in [self.rheology, self.radiogenics, self.cooling_model]:
            model.clear_state()

    def set_state(self, temperature: FloatArray = None, pressure: FloatArray = None,
                  force_update_strength: bool = True):
        """ Set the layer's temperature and update all related properties.

        Parameters
        ----------
        temperature : FloatArray = None
            Temperature of layer [K]
        pressure : FloatArray = None
            Pressure of layer [Pa]
        force_update_strength : bool = True
            If True, <layer>.update_strength() will be called

        """

        if temperature is not None:
            new_shape, temperature = reshape_help(temperature, self.world.global_shape,
                                                  call_locale=f'{self}.set_state.temperature')
            if new_shape:
                self.world.change_shape(new_shape)

            self._temperature = temperature

        if pressure is not None:
            new_shape, pressure = reshape_help(pressure, self.world.global_shape,
                                               call_locale=f'{self}.set_state.pressure')
            if new_shape:
                self.world.change_shape(new_shape)
            self._pressure = pressure

        if force_update_strength:
            # Temperature and pressure will change the strength of the layer and all of its dependencies
            self.update_strength()

    def set_temperature(self, temperature: FloatArray):
        """ Wrapper for PhysicsLayer.set_state

        Parameters
        ----------
        temperature : FloatArray
            Temperature of layer [K]

        See Also
        --------
        TidalPy.structures.layers.layers.PhysicsLayer.set_state
        """

        self.set_state(temperature=temperature, pressure=None)

    def set_pressure(self, pressure: FloatArray):
        """ Wrapper for PhysicsLayer.set_state

        Parameters
        ----------
        pressure : FloatArray
            Pressure of layer [Pa]

        See Also
        --------
        TidalPy.structures.layers.layers.PhysicsLayer.set_state
        """

        self.set_state(temperature=None, pressure=pressure)

    def set_strength(self, viscosity: FloatArray = None, shear_modulus: FloatArray = None):
        """ Manual set the viscosity and shear modulus of the layer, independent of temperature.

        This method by-passes the self.viscosity_func and allows the user to manually set the viscosity and shear
        of the layer.

        Future : TODO
        ------
        Currently it does not change the layer's temperature. In the future an effective temperature can be
            (optionally) calculated. Make sure to implement this at the Rheology class level.

        Parameters
        ----------
        viscosity : FloatArray
            The new viscosity of the layer in [Pa s]
        shear_modulus : FloatArray
            The new shear modulus of the layer in [Pa]
        """

        # Pass the new strength to the rheology class
        self.rheology.set_state(viscosity, shear_modulus, called_by_layer=True)

        self.update_strength(already_called_rheology=True)

    def update_radiogenics(self, force_calculation: bool = False) -> np.ndarray:
        """ Calculate the radiogenic heating rate within this layer, requires the world.time to be set.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.

        Returns
        -------
        radiogenic_heating : np.ndarray
            Radiogenic heating rate in [Watts]
        """

        try:
            self.radiogenics.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error

        # New radiogenic heating will change the temperature derivative
        self.calc_temperature_derivative(force_calculation=False)

        return self.radiogenic_heating

    def update_cooling(self, force_calculation: bool = False, call_deriv_calc: bool = True):
        """ Calculate parameters related to cooling of the layer, be it convection or conduction

        Using the self.cooling class and the current temperature and viscosity, this method will calculate convection
            parameters. Some of these may be zeros if convection is forced off in the planet's configuration.

        By default this method will ignore parameter missing errors to avoid initialization errors. This behaviour
            can be altered by the force_calculation flag.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.
        call_deriv_calc : bool = True
            If True, then this method will attempt to call self.calc_temperature_derivative
        """

        try:
            self.cooling_model.calculate()
        except ParameterMissingError as error:
            if force_calculation:
                raise error

        if self.is_top_layer and self.world.orbit is not None:
            # TODO: Add some sort of convergence here?
            # Tell the planet that there is a new surface flux and thus a new surface temperature
            self.world.update_surface_temperature()

        if call_deriv_calc:
            # New cooling will change the temperature derivative
            self.calc_temperature_derivative(force_calculation=False)

    def update_tides(self, force_calculation: bool = False, call_deriv_calc: bool = True):
        """ Calculate, and update the layer's state variables with tidal heating and torques based on its current state

        Using the self.tides class and the layer's current temperature, viscosity, shear modulus, and frequencies;
            the layer will attempt to calculate a new effective rigidity, complex compliance, and complex love number,
            ultimately it will find tidal heating and torque then inform the planet that it should update its own
            global calculations.

        Viscosity, Shear Modulus, and Orbital Freq are required inputs to most rheological models. If the planet is not
            forced into spin-sync, then Spin Freq is also required. If a missing parameter is raised then the method
            will exit unless force_calculation is set to True.

        Parameters
        ----------
        force_calculation : bool = False
            Flag that will direct the method to re-raise any missing parameter errors that may otherwise be ignored.
        call_deriv_calc : bool = True
            If True, then this method will attempt to call self.calc_temperature_derivative

        """

        if self.is_tidal:
            try:
                self.world.update_tides()
            except ParameterMissingError as error:
                if force_calculation:
                    raise error

        if call_deriv_calc:
            self.calc_temperature_derivative(force_calculation=False)

    def update_strength(self, already_called_rheology: bool = False):
        """ Calculate viscosity and shear modulus and update the layer's state with new values.

        Rheology class handles the change in pre-melt shear modulus and viscosity as well as the impact of partial
            melting. It will then calculate the complex compliances which will impact the global tides.

        The pressure of the layer will be used, if present, in the calculations (assuming that the viscosity and shear
            functions actually depend upon pressure). Otherwise, a pressure of zero will be assumed.

        Parameters
        ----------
        already_called_rheology : bool = False
            Some other methods of this class may call this method and may have already performed some of the tasks
                that this method would normally perform.

        """

        if not already_called_rheology:
            # Tell rheology class to update strength (calculation of effective viscosity, shear, and complex compliance)
            self.rheology.update_strength(called_from_layer=True)

        # Changing the strength will change the cooling properties. Don't have it call derivative calculation yet.
        self.update_cooling(call_deriv_calc=False)

        # Changing the strength will change tides
        self.update_tides(call_deriv_calc=False)

        # The above changes will impact the temperature derivative calculation. So let's call that now.
        self.calc_temperature_derivative(force_calculation=False)

    def calc_temperature_derivative(self, force_calculation: bool = False) -> np.ndarray:
        """ Calculate the change in temperature within this layer over time.

        dT/dt = (heating_in - cooling_out) / (M * c_p * (St + 1) * Tr)

        Returns
        -------
        dT/dt : np.ndarray
            The derivative of temperature with respect to time.
        """

        # Initialize total heating based on if there is a layer warming this one from below
        if self.layer_below is None:
            total_heating = np.zeros_like(self.temperature)
        elif self.layer_below.heat_flux is None:
            total_heating = np.zeros_like(self.temperature)
        else:
            total_heating = self.layer_below.heat_flux * self.surface_area_inner

        try:
            # Heat sources are added in the reinit() method.
            for heat_source_func in self.heat_sources:
                heat_source = heat_source_func()
                if heat_source is None:
                    raise ParameterMissingError
                total_heating += heat_source

            # Cooling is calculated by the cooling_model class
            total_cooling = self.cooling

            self._deriv_temperature = (total_heating - total_cooling) / \
                                     (self.mass * self.specific_heat * (self.stefan + 1.) * self.temp_ratio)
        except ParameterMissingError as error:
            if force_calculation:
                raise error

        return self.deriv_temperature

    def geotherm(self, avg_temperature: float = None):
        """ Calculates layer's geotherm based on an average temperature (or layer's current temperature)

        Returns
        -------
        temperature_profile : np.ndarray
            numpy array of the adiabatic temperature profile
        """

        if avg_temperature is None:
            avg_temperature = self.temperature

        temperature_profile = avg_temperature

        return temperature_profile

    # State properties
    @property
    def deriv_temperature(self) -> np.ndarray:
        """ Time Derivative of Temperature [T s-1] """
        return self._deriv_temperature

    @deriv_temperature.setter
    def deriv_temperature(self, value):
        raise ImproperPropertyHandling


    # Inner-scope properties
    # # Rheology Class
    @property
    def viscosity(self):
        """ Solid Viscosity of the Layer [Pa s]

        Wrapper for `<layer>.rheology.postmelt_viscosity` """

        return self.rheology.postmelt_viscosity

    @viscosity.setter
    def viscosity(self, value):
        self.rheology.postmelt_viscosity = value

    @property
    def liquid_viscosity(self):
        """ Liquid Viscosity of the Layer [Pa s]

        Wrapper for `<layer>.rheology.liquid_viscosity` """

        return self.rheology.liquid_viscosity

    @liquid_viscosity.setter
    def liquid_viscosity(self, value):
        self.rheology.liquid_viscosity = value

    @property
    def shear_modulus(self):
        """ Shear Modulus (Shear Rigidity) of the Layer [Pa]

        Wrapper for `<layer>.rheology.postmelt_shear_modulus` """

        return self.rheology.postmelt_shear_modulus

    @shear_modulus.setter
    def shear_modulus(self, value):
        self.rheology.postmelt_shear_modulus = value

    @property
    def compliance(self):
        """ Shear Compliance (Inverse of Shear Rigidity) of the Layer [Pa-1]

        Wrapper for `<layer>.rheology.postmelt_compliance` """

        return self.rheology.postmelt_compliance

    @compliance.setter
    def compliance(self, value):
        self.rheology.postmelt_compliance = value

    @property
    def complex_compliances(self):
        """ Complex Shear Compliance of the Layer [Pa-1]

        This is a complex number which includes information on the layer's ability to dissipate shear energy.

        Wrapper for `<layer>.rheology.complex_compliances` """

        return self.rheology.complex_compliances

    @complex_compliances.setter
    def complex_compliances(self, value):
        self.rheology.complex_compliances = value

    @property
    def melt_fraction(self):
        """ Melt Fraction of the Layer

        Wrapper for `<layer>.rheology.melt_fraction` """

        return self.rheology.melt_fraction

    @melt_fraction.setter
    def melt_fraction(self, value):
        self.rheology.melt_fraction = value

    # # Radiogenics Class
    @property
    def radiogenic_heating(self):
        """ Radiogenic Heating Rate [W]

        Wrapper for `<layer>.radiogenics.heating`

        Notes:
        Calculated at the <world>.time (or, equivalently, <layer>.time)
        """

        return self.radiogenics.heating

    @radiogenic_heating.setter
    def radiogenic_heating(self, value):
        self.radiogenics.heating = value

    # # Cooling Class
    @property
    def cooling(self):
        """ Layer Cooling Rate [W]

        Wrapper for `<layer>.cooling_model.cooling`

        Notes
        -----
        Cooling module determines if the layer is convecting or conduction (depending on user-provided configurations).
        """

        return self.cooling_model.cooling

    @cooling.setter
    def cooling(self, value):
        self.cooling_model.cooling = value

    @property
    def cooling_flux(self):
        """ Layer Cooling Flux [W m-2]

        Wrapper for `<layer>.cooling_model.cooling_flux`

        Notes
        -----
        Cooling module determines if the layer is convecting or conduction (depending on user-provided configurations).
        """

        return self.cooling_model.cooling_flux

    @cooling_flux.setter
    def cooling_flux(self, value):
        self.cooling_model.cooling_flux = value

    @property
    def boundary_layer_thickness(self):
        """ Layer's Thermal Boundary Layer Thickness [m]

        Wrapper for `<layer>.cooling_model.boundary_layer_thickness`

        Notes
        -----
        Cooling module determines if the layer is convecting or conduction (depending on user-provided configurations).
        If the layer is set to only conduct, then the boundary layer thickness = 0.5 layer thickness.
        """

        return self.cooling_model.boundary_layer_thickness

    @boundary_layer_thickness.setter
    def boundary_layer_thickness(self, value):
        self.cooling_model.boundary_layer_thickness = value

    @property
    def rayleigh(self):
        """ Layer's Thermal Rayleigh Number

        Wrapper for `<layer>.cooling_model.rayleigh`

        Notes
        -----
        Cooling module determines if the layer is convecting or conduction (depending on user-provided configurations).
        If the layer is not convecting then the Rayleigh number will be 0
        """

        return self.cooling_model.rayleigh

    @rayleigh.setter
    def rayleigh(self, value):
        self.cooling_model.rayleigh = value

    @property
    def nusselt(self) -> np.ndarray:
        """ Layer's Thermal Nusselt Number

        Wrapper for `<layer>.cooling_model.nusselt`

        Notes
        -----
        Cooling module determines if the layer is convecting or conduction (depending on user-provided configurations).
        If the layer is not convecting then the Nusselt number will be 1
        """

        return self.cooling_model.nusselt

    @nusselt.setter
    def nusselt(self, value):
        self.cooling_model.nusselt = value


    # Outer-scope properties
    # # Tides Class
    @property
    def tidal_heating(self):
        """ Layer's Tidal Heating

        Wrapper for `<layer>.<world>.tides.tidal_heating_by_layer[self]`
        """

        return self.world.tides.tidal_heating_by_layer[self]

    @tidal_heating.setter
    def tidal_heating(self, value):
        raise OuterscopePropertySetError


    # Aliased properties
    @property
    def blt(self):
        """ Alias for PhysicsLayer.boundary_layer_thickness """
        return self.boundary_layer_thickness

    @blt.setter
    def blt(self, value):
        self.boundary_layer_thickness = value