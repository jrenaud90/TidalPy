from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

from TidalPy.exceptions import (AttributeNotSetError, ConfigPropertyChangeError, ImproperPropertyHandling,
                                MissingAttributeError, OuterscopePropertySetError)

from .basic import LayerBase

from TidalPy.logger import get_logger
log = get_logger("TidalPy")


if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray, NoneType

    from ..world_types import TidalWorldType


class PhysicsLayer(LayerBase):
    """ PhysicsLayer
    Layer object to store parameters geometric and physical properties calculated by TidalPy based on a user-provided
        configuration dictionary. PhysicsLayer contains additional properties and methods used to perform various
        thermal and tidal calculations.

    Notes:
    .. Should not be used with Burnman calculated layers (see BurnmanLayer instead)

    See Also
    --------
    TidalPy.structures.layers.LayerBase
    TidalPy.structures.layers.BurnmanLayer
    """

    layer_class = 'physics'

    def __init__(
        self, layer_name: str, layer_index: int, world: 'TidalWorldType', layer_config: dict,
        is_top_layer: bool, initialize: bool = True
        ):
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
        self._temperature_time_derivative = None

        # Configuration properties
        self._use_pressure_in_strength_calc = False
        # Model holders - reinit in reinit()
        self._rheology = None       # type: Union['Rheology', NoneType]
        self._cooling_model = None  # type: Union['CoolingModel', NoneType]
        self._radiogenics = None    # type: Union['Radiogenics', NoneType]

        self.heat_sources = list()

        # Material Properties
        # TODO: Should all of these be private properties with a dedicated setter that makes a call to update_thermal()?
        # These are calculated by BurnMan in the BurnmanLayer version, can be set by user/config here.
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

    def reinit(self, initial_init: bool = False, initialize_geometry: bool = True):
        """ Reinitialize the physical object by pulling in any potentially new configurations

        Parameters
        ----------
        initial_init : bool = False
            Set to `True` for the first time an instance is created.
        initialize_geometry : bool = False
            Set to `True` if the set_geometry method should be called from within reinit
        """

        # Base class's reinit is called *after* the geometry is set (so that tidal volume fraction is set correctly)
        super().reinit(initial_init=initial_init, initialize_geometry=initialize_geometry)

        # Material properties that might have been affected by new configuration files
        self.static_shear_modulus = self.config['shear_modulus']
        self.heat_fusion = self.config['heat_fusion']
        self.thermal_conductivity = self.config['thermal_conductivity']
        self.thermal_diffusivity = self.config['thermal_diffusivity']
        self.thermal_expansion = self.config['thermal_expansion']
        self.temp_ratio = self.config['boundary_temperature_ratio']
        self.stefan = self.config['stefan']

        # Setup Models
        from TidalPy.rheology import Rheology
        from TidalPy.radiogenics.radiogenics import Radiogenics
        from TidalPy.cooling import CoolingModel

        self._rheology = Rheology(self, store_config_in_layer=True)
        self._radiogenics = Radiogenics(self, store_config_in_layer=True)
        self._cooling_model = CoolingModel(self, store_config_in_layer=True)

        # Find all heating sources
        self.heat_sources = list()
        if 'off' not in self.radiogenics.model:
            self.heat_sources.append(lambda: self.radiogenic_heating)
        if self.is_tidal:
            self.heat_sources.append(lambda: self.tidal_heating)

    def update_cooling(self, force_update: bool = False):
        """ Calculate parameters related to cooling of the layer, be it convection or conduction

        Using the self.cooling class and the current temperature and viscosity, this method will calculate convection
            parameters. Some of these may be zeros if convection is forced off in the planet's configuration.

        By default this method will ignore parameter missing errors to avoid initialization errors. This behaviour
            can be altered by the force_calculation flag.

        Parameters:
        force_update : bool = False
            If `False`, then any (expected) errors are raised during the cooling calculation will be ignored.

        """

        # Update the layer's cooling
        cooling_calculated = False
        try:
            self.cooling_model.calculate()
            cooling_calculated = True
        except MissingAttributeError as e:
            if force_update:
                raise e

        # If the layer's cooling changed, update the thermal equilibrium
        if self._cooling_model.model.lower() != 'off' and cooling_calculated:
            self.internal_thermal_equilibrium_changed(called_from_cooling=True)

    def time_changed(self):
        """ The time has changed. Make any necessary updates. """

        super().time_changed()

        self.radiogenics.calculate()

        if self.radiogenics.model.lower() != 'off':
            # Tell the layer that its internal thermal equilibrium has changed.
            self.internal_thermal_equilibrium_changed()

    def temperature_pressure_changed(self):
        """ The temperature and/or pressure of the layer has changed. Make any necessary updates. """

        super().temperature_pressure_changed()

        # Tell the rheology class that the temperature and/or pressure has changed.
        self.rheology.temperature_pressure_changed()

    def strength_changed(self):
        """ The viscosity and/or shear modulus of the layer has changed. Make any necessary updates. """

        super().strength_changed()

        # Tell the rheology class that the viscosity and/or rigidity has changed.
        self.rheology.strength_changed()

    def internal_thermal_equilibrium_changed(self, called_from_cooling: bool = False):
        """ The internal heating / cooling of the layer has changed. Make any necessary updates.

        Parameters
        ----------
        called_from_cooling : bool = False
            Flag to avoid recursive loops between surface temperature and cooling.
        """

        super().internal_thermal_equilibrium_changed()

        if self.is_top_layer:
            # Tell the world that the surface temperature needs to update
            self.world.update_surface_temperature(called_from_cooling=called_from_cooling)

        if self.cooling is not None:
            # Calculate the temperature derivative of the layer
            self.calc_temperature_derivative()

    def surface_temperature_changed(self, called_from_cooling: bool = False):
        """ Surface temperature has changed - Perform any calculations that may have also changed.

        Parameters
        ----------
        called_from_cooling : bool = False
            Flag to avoid recursive loops between surface temperature and cooling.
        """

        super().surface_temperature_changed()

        # Update cooling model based on new surface temperature.
        if not called_from_cooling:
            self.update_cooling()

    def tidal_frequencies_changed(self, collapse_tidal_modes: bool = True):
        """ The tidal frequencies have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        super().tidal_frequencies_changed(collapse_tidal_modes=collapse_tidal_modes)

        # Tell the rheology to update its complex compliances if applicable.
        self.rheology.tidal_frequencies_changed(collapse_tidal_modes=collapse_tidal_modes)

    def complex_compliances_changed(self, collapse_tidal_modes: bool = True):
        """ The complex compliances have changed. Make any necessary updates.

        Parameters
        ----------
        collapse_tidal_modes : bool = True
            If `True`, then the world will tell its tides model to collapse tidal modes.
        """

        # This is called from bottom-to-top starting in the ComplexCompliances class inside Rheology.
        self.world.complex_compliances_changed(collapse_tidal_modes=collapse_tidal_modes)

    def clear_state(self, clear_pressure: bool = False):

        super().clear_state(clear_pressure=clear_pressure)

        # State Derivatives
        self._temperature_time_derivative = None

        # Clear the state of all inner-scope methods
        for model in [self.rheology, self.radiogenics, self.cooling_model]:
            model.clear_state()

    def set_strength(self, viscosity: 'FloatArray' = None, shear_modulus: 'FloatArray' = None):
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
        self.rheology.set_state(viscosity, shear_modulus)

    def set_state(
        self, temperature: 'FloatArray' = None, pressure: 'FloatArray' = None, viscosity: 'FloatArray' = None,
        shear_modulus: 'FloatArray' = None
        ):
        """ Set the layer's state properties

        Parameters
        ----------
        temperature : FloatArray = None
            New dynamic temperature for the layer [K].
        pressure : FloatArray = None
            New dynamic pressure for the layer [Pa].
        viscosity : FloatArray = None
            The new viscosity of the layer in [Pa s]
        shear_modulus : FloatArray = None
            The new shear modulus of the layer in [Pa]

        """

        # Check if too much information was set.
        temperature_change = temperature is not None
        pressure_change = pressure is not None
        viscosity_change = viscosity is not None
        shear_modulus_change = shear_modulus is not None

        if (temperature_change or pressure_change) and (viscosity_change or shear_modulus_change):
            log.warning(
                f'Trying to set TP and strength at the same time for layer {self}. '
                f'Only setting the strength changes.'
                )
            temperature_change = False
            pressure_change = False

        # Update the temperature and pressure state using the parent classes method.
        if temperature_change or pressure_change:
            super().set_state(temperature=temperature, pressure=pressure)

        # If the strength is being changed, use this layer classe's set_strength method.
        if viscosity_change or shear_modulus_change:
            self.set_strength(viscosity=viscosity, shear_modulus=shear_modulus)

    def calc_temperature_derivative(self, force_calculation: bool = False) -> FloatArray:
        """ Calculate the change in temperature within this layer over time.

        .. math:: dT / dt = (Q_{ \text{In} } - Q_{ \text{Out} }) / (M * c_{ \text{p} } * (\text{St} + 1)  * T_{r}

        Parameters
        ----------
        force_calculation : bool = False
            If `True`, then the method will raise any errors encountered. Otherwise method will ignore errors
                but the temperature_time_derivative will likely not be set.

        Returns
        -------
        dT/dt : FloatArray
            The derivative of temperature with respect to time [K s-1].
        """

        # Initialize the total heating rate of this layer based on if there is a layer warming this one from below.
        if self.layer_below is None:
            total_heating = np.zeros_like(self.temperature)
        elif self.layer_below.cooling_flux is None:
            total_heating = np.zeros_like(self.temperature)
        else:
            total_heating = self.layer_below.cooling_flux * self.surface_area_inner

        # Add additional heat sources that are setup in this layer.
        try:
            # Heat sources are added in the reinit() method.
            for heat_source_func in self.heat_sources:
                heat_source = heat_source_func()
                if heat_source is None:
                    log.warning(f'One or more heat sources were not set for {self}.')
                else:
                    total_heating += heat_source

            # Cooling is calculated by the cooling_model class
            if self.cooling is None:
                log.warning(f'Cooling was not set for layer {self}.')
                self._temperature_time_derivative = None
            else:
                total_cooling = self.cooling
                self._temperature_time_derivative = \
                    (total_heating - total_cooling) / \
                    (self.mass * self.specific_heat * (self.stefan + 1.) * self.temp_ratio)

        except AttributeNotSetError as error:
            self._temperature_time_derivative = None
            if force_calculation:
                raise error
        except ValueError as error:
            self._temperature_time_derivative = None
            if force_calculation:
                raise error
        except TypeError as error:
            self._temperature_time_derivative = None
            if force_calculation:
                raise error

        return self.temperature_time_derivative

    def geotherm(self, avg_temperature: float = None):
        """ Calculates layer's geotherm based on an average temperature (or layer's current temperature)

        Returns
        -------
        temperature_profile : np.ndarray
            numpy array of the adiabatic temperature profile
        """

        if avg_temperature is None:
            if type(self.temperature) == np.ndarray:
                avg_temperature = np.mean(self.temperature)
            else:
                avg_temperature = self.temperature

        temperature_profile = avg_temperature

        return temperature_profile

    # # State properties
    @property
    def temperature_time_derivative(self) -> FloatArray:
        """ Time Derivative of Temperature [K s-1] """
        return self._temperature_time_derivative

    @temperature_time_derivative.setter
    def temperature_time_derivative(self, value):
        raise ImproperPropertyHandling

    # # Configuration properties
    @property
    def use_pressure_in_strength_calc(self) -> bool:
        """ Flag for if pressure is used in viscosity and rigidity calculations """
        return self._use_pressure_in_strength_calc

    @use_pressure_in_strength_calc.setter
    def use_pressure_in_strength_calc(self, value):
        raise ConfigPropertyChangeError

    # Model storage
    @property
    def rheology(self) -> 'Rheology':
        """ Rheology class instance, used to calculate viscosity, rigidity, partial melting, and complex compliance """
        return self._rheology

    @rheology.setter
    def rheology(self, value):
        raise ConfigPropertyChangeError

    @property
    def cooling_model(self) -> 'CoolingModel':
        """ Cooling class instance, used to calculate cooling rates """
        return self._cooling_model

    @cooling_model.setter
    def cooling_model(self, new_cooling_model: 'CoolingModel'):
        raise ConfigPropertyChangeError

    @property
    def radiogenics(self) -> 'Radiogenics':
        """ Radiogenics class instance, used to calculate radiogenic heating based on the current time """
        return self._radiogenics

    @radiogenics.setter
    def radiogenics(self, new_radiogenics: 'Radiogenics'):
        raise ConfigPropertyChangeError

    # Inner-scope properties
    # # Rheology Class
    @property
    def viscosity(self):
        """ Solid Viscosity of the Layer [Pa s]

        Wrapper for `<layer>.rheology.postmelt_viscosity` """

        return self.rheology.postmelt_viscosity

    @viscosity.setter
    def viscosity(self, new_viscosity: FloatArray):
        self.set_strength(viscosity=new_viscosity)

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
    def shear_modulus(self, new_shear_modulus: FloatArray):
        self.set_strength(shear_modulus=new_shear_modulus)

    @property
    def compliance(self):
        """ Shear Compliance (Inverse of Shear Rigidity) of the Layer [Pa-1]

        Wrapper for `<layer>.rheology.postmelt_compliance` """

        return self.rheology.postmelt_compliance

    @compliance.setter
    def compliance(self, new_compliance: FloatArray):
        self.set_strength(shear_modulus=new_compliance**(-1))

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

        if self.world.tides is None:
            return None
        else:
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
