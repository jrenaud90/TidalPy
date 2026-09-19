# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's gas layer class.

GasLayer extends PhysicsLayer with ideal-gas thermodynamics (adiabatic lapse rate, scale height, ideal-gas
pressure, sound speed). It has no phase changes and no cooling or radiogenics sub-models.
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr, make_unique

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer, c_tidal_scale_method_from_name
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# GasLayer
# =====================================================================================================================
cdef class GasLayer(PhysicsLayer):
    """Ideal-gas layer: PhysicsLayer plus adiabatic lapse rate, scale height, ideal-gas pressure, and sound speed.

    No phase changes, cooling, or radiogenics sub-models are available (use SolidLiquidLayer for those).

    Parameters
    ----------
    name : str
        Layer name.
    layer_index : int
        Zero-based position; innermost layer = 0.
    radius_inner : float
        Inner boundary radius [m].
    radius_outer : float
        Outer boundary radius [m].
    mass : float
        Total layer mass [kg].
    material_name : str, optional
        Material identifier. Default ``""``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        Dimensionless tidal heating scale. Default ``1.0``.
    shear_modulus_static : float, optional
        Unrelaxed shear modulus [Pa]. Default ``0.0``.
    bulk_modulus_static : float, optional
        Unrelaxed bulk modulus [Pa]. Default ``0.0``.
    shear_viscosity_static : float, optional
        Reference dynamic shear viscosity [Pa·s]. Default NaN (unset).
    bulk_viscosity_static : float, optional
        Reference dynamic bulk viscosity [Pa·s]. Default NaN (unset).
    love_number_k : complex, optional
        Potential Love number k (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number h (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number l (placeholder). Default ``0+0j``.
    mean_molecular_weight : float, optional
        Mean molecular weight of the gas [kg/mol]. Default ``2e-3`` (H₂).
    adiabatic_index : float, optional
        Ratio of specific heats γ = c_p/c_v [dimensionless]. Default ``1.4``.
    reference_temperature : float, optional
        Reference temperature [K]. Default ``300.0``.
    reference_density : float, optional
        Reference density [kg/m³]. Default ``1.0``.
    tidal_scale_method : str, optional
        How the layer's share of the world's tidal heating is set. Default ``"user_provided"``.
    is_solid : bool, optional
        True marks the layer solid for the radial Love-number solver. Default ``False``: a gas carries no
        shear stress, so it is solved as a liquid.
    is_static : bool, optional
        Use the static (no inertia) approximation in the radial solver. Default ``True``.
    is_incompressible : bool, optional
        Use the incompressible approximation in the radial solver. Default ``False``.
    temperature : float, optional
        Layer temperature [K] at which its viscosity and melt models are evaluated. Default ``0.0``, the cold
        rigid limit of the viscosity laws.
    shear_modulus_pressure_derivative : float, optional
        Pressure derivative of the static shear modulus [Pa/Pa]. Default ``0.0``.
    shear_modulus_temperature_derivative : float, optional
        Temperature derivative of the static shear modulus [Pa/K]. Default ``0.0``.
    shear_modulus_reference_temperature : float, optional
        Temperature [K] at which ``shear_modulus_static`` applies. ``None`` keeps the default of 300 K.
    use_thermal_eos : bool, optional
        Pass the temperature to the EOS model, so the density and bulk modulus depend on it. Default ``False``.

    Assumptions
    -----------
    - Spherically symmetric layer geometry.
    - The gas is ideal; real-gas corrections are not included.
    - The universal gas constant comes from the TidalPy global configuration.
    """

    def __cinit__(self, *args, **kwargs):
        self._gas_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner,
            double radius_outer,
            double mass,
            str    material_name        = "",
            cpp_bool is_tidal           = True,
            double tidal_scale          = 1.0,
            double shear_modulus_static = 0.0,
            double bulk_modulus_static  = 0.0,
            double shear_viscosity_static = d_NAN,
            double bulk_viscosity_static = d_NAN,
            complex love_number_k        = 0+0j,
            complex love_number_h        = 0+0j,
            complex love_number_l        = 0+0j,
            double mean_molecular_weight = 2.0e-3,
            double adiabatic_index       = 1.4,
            double reference_temperature = 300.0,
            double reference_density  = 1.0,
            str    tidal_scale_method = "user_provided",
            cpp_bool is_solid          = False,
            cpp_bool is_static         = True,
            cpp_bool is_incompressible = False,
            double temperature           = 0.0,
            double shear_modulus_pressure_derivative    = 0.0,
            double shear_modulus_temperature_derivative = 0.0,
            shear_modulus_reference_temperature         = None,
            cpp_bool use_thermal_eos     = False):
        cdef c_GasConfig config
        config.name                 = name.encode("utf-8")
        config.layer_index          = layer_index
        config.radius_inner         = radius_inner
        config.radius_outer         = radius_outer
        config.mass                 = mass
        config.material_name        = material_name.encode("utf-8")
        config.is_tidal             = is_tidal
        config.tidal_scale          = tidal_scale
        config.tidal_scale_method   = c_tidal_scale_method_from_name(tidal_scale_method.encode("utf-8"))
        config.shear_modulus_static = shear_modulus_static
        config.bulk_modulus_static  = bulk_modulus_static
        config.shear_viscosity_static = shear_viscosity_static
        config.bulk_viscosity_static  = bulk_viscosity_static
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.is_solid              = is_solid
        config.is_static             = is_static
        config.is_incompressible     = is_incompressible
        config.temperature       = temperature
        config.shear_modulus_pressure_derivative    = shear_modulus_pressure_derivative
        config.shear_modulus_temperature_derivative = shear_modulus_temperature_derivative
        # None keeps the C++ default reference temperature.
        if shear_modulus_reference_temperature is not None:
            config.shear_modulus_reference_temperature = <double>shear_modulus_reference_temperature
        config.use_thermal_eos   = use_thermal_eos
        config.mean_molecular_weight = mean_molecular_weight
        config.adiabatic_index       = adiabatic_index
        config.reference_temperature = reference_temperature
        config.reference_density     = reference_density
        # make_unique owns the allocation; ownership then moves into the base-typed member
        # (Cython cannot assign a unique_ptr[Derived] to a unique_ptr[Base] directly).
        cdef unique_ptr[c_GasLayer] built = make_unique[c_GasLayer](config)
        self._gas_ptr     = built.get()
        self._physics_ptr = <c_PhysicsLayer*>self._gas_ptr
        self._layer_ptr.reset(<c_BaseLayer*>built.release())
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()

    def __dealloc__(self):
        self._gas_ptr     = NULL  # base's unique_ptr owns the C++ object
        self._physics_ptr = NULL

    @staticmethod
    cdef GasLayer _view(c_GasLayer* ptr, object world):
        cdef GasLayer v = GasLayer.__new__(GasLayer)
        v._gas_ptr      = ptr
        v._physics_ptr  = <c_PhysicsLayer*>ptr
        v._init_view(<c_BaseLayer*>ptr, world)
        return v

    # ------------------------------------------------------------------------------------------------------------------
    # Gas properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def mean_molecular_weight(self) -> float:
        """Mean molecular weight of the gas [kg/mol]."""
        return self._gas_ptr.get_mean_molecular_weight()

    @property
    def adiabatic_index(self) -> float:
        """Ratio of specific heats γ = c_p/c_v [dimensionless]."""
        return self._gas_ptr.get_adiabatic_index()

    @property
    def reference_temperature(self) -> float:
        """Reference temperature [K]."""
        return self._gas_ptr.get_reference_temperature()

    @property
    def reference_density(self) -> float:
        """Reference density [kg/m³]."""
        return self._gas_ptr.get_reference_density()

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_adiabatic_lapse_rate(self, double gravity) -> float:
        """Dry adiabatic lapse rate [K/m] = g * (γ-1) * M / (γ * R); 0.0 for invalid inputs.

        Parameters
        ----------
        gravity : float
            Gravitational acceleration [m/s²].

        Returns
        -------
        float
            Adiabatic lapse rate [K/m].
        """
        return self._gas_ptr.calc_adiabatic_lapse_rate(gravity)

    def calc_scale_height(self, double temperature,
                          double gravity) -> float:
        """Barometric scale height [m] = R * T / (g * M); 0.0 for invalid inputs.

        Parameters
        ----------
        temperature : float
            Temperature [K].
        gravity : float
            Gravitational acceleration [m/s²].

        Returns
        -------
        float
            Scale height [m].
        """
        return self._gas_ptr.calc_scale_height(temperature, gravity)

    def calc_pressure_ideal_gas(self, double temperature,
                                double density) -> float:
        """Ideal-gas pressure [Pa] = ρ * R * T / M; 0.0 for invalid inputs.

        Parameters
        ----------
        temperature : float
            Temperature [K].
        density : float
            Gas density [kg/m³].

        Returns
        -------
        float
            Pressure [Pa].
        """
        return self._gas_ptr.calc_pressure_ideal_gas(temperature, density)

    def calc_sound_speed(self, double temperature) -> float:
        """Adiabatic sound speed [m/s] = sqrt(γ * R * T / M); 0.0 for invalid inputs.

        Parameters
        ----------
        temperature : float
            Temperature [K].

        Returns
        -------
        float
            Sound speed [m/s].
        """
        return self._gas_ptr.calc_sound_speed(temperature)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS): the PhysicsLayer keys plus the gas parameters."""
        d = PhysicsLayer.get_config_dict(self)
        d["mean_molecular_weight_kg_mol"] = self._gas_ptr.get_mean_molecular_weight()
        d["adiabatic_index"]              = self._gas_ptr.get_adiabatic_index()
        d["reference_temperature_k"]      = self._gas_ptr.get_reference_temperature()
        d["reference_density_kg_m3"]      = self._gas_ptr.get_reference_density()
        return d
