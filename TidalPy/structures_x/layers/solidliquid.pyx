# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
solidliquid.pyx
Cython/Python wrapper for TidalPy's solid/liquid layer class.

SolidLiquidLayer: extends PhysicsLayer with thermal properties, Arrhenius
viscosity, melt-fraction tracking, and optional cooling/radiogenics sub-models.
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer, c_tidal_scale_method_from_name
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers
from TidalPy.cooling_x.cooling cimport CoolingBase
from TidalPy.radiogenics_x.radiogenics cimport RadiogenicsBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# SolidLiquidLayer
# =====================================================================================================================
cdef class SolidLiquidLayer(PhysicsLayer):
    """Thermo-mechanical layer with phase-change tracking, Arrhenius viscosity,
    and optional cooling/radiogenics sub-models.

    Extends PhysicsLayer with:

    - Temperature-dependent melt fraction (linear between solidus and liquidus).
    - Arrhenius viscosity with pressure and partial-melt corrections.
    - Melt-fraction-reduced shear modulus.
    - Thermal conductivity, diffusivity, and adiabatic temperature gradient.
    - Conductive heat flux through the layer.
    - Radiogenic heating via an optional RadiogenicsBase sub-model (Phase 7).
    - Convective/conductive cooling via an optional CoolingBase sub-model (Phase 6).

    Parameters
    ----------
    name : str
        Human-readable layer name.
    layer_index : int
        Zero-based position; innermost layer = 0.
    radius_inner_m : float
        Inner boundary radius [m].
    radius_outer_m : float
        Outer boundary radius [m].
    mass_kg : float
        Total layer mass [kg].
    material_name : str, optional
        Material identifier. Default ``""``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        Dimensionless tidal heating scale. Default ``1.0``.
    shear_modulus_static_pa : float, optional
        Unrelaxed (zero-temperature) shear modulus [Pa]. Default ``0.0``.
    bulk_modulus_static_pa : float, optional
        Unrelaxed bulk modulus [Pa]. Default ``0.0``.
    shear_viscosity_static_pas : float, optional
        Reference shear viscosity at reference_temperature_k and P=0 [Pa·s]. Default ``0.0``.
    bulk_viscosity_static_pas : float, optional
        Reference bulk viscosity [Pa·s]. Default ``0.0``.
    love_number_k : complex, optional
        Potential Love number k (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number h (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number l (placeholder). Default ``0+0j``.
    thermal_conductivity_ref_w_mk : float, optional
        Reference thermal conductivity [W/(m·K)]. Default ``4.0``.
    thermal_expansion_ref_1_k : float, optional
        Reference thermal expansion coefficient [1/K]. Default ``3e-5``.
    heat_capacity_ref_j_kgk : float, optional
        Reference specific heat capacity [J/(kg·K)]. Default ``1200.0``.
    activation_energy_j_mol : float, optional
        Arrhenius activation energy [J/mol]. Default ``300e3``.
    activation_volume_m3_mol : float, optional
        Arrhenius activation volume [m³/mol]. Default ``5e-6``.
    solidus_temperature_k : float, optional
        Solidus temperature [K]. Default ``1600.0``.
    liquidus_temperature_k : float, optional
        Liquidus temperature [K]. Default ``2000.0``.
    melt_fraction_exponent : float, optional
        Exponent in melt-fraction parameterization. Default ``1.0``.
    reference_density_kg_m3 : float, optional
        Reference density for thermal diffusivity [kg/m³]. Default ``3500.0``.
    reference_temperature_k : float, optional
        Reference temperature for Arrhenius viscosity [K]. Default ``1600.0``.
    melt_viscosity_reduction : float, optional
        Exponential melt-viscosity reduction coefficient. Default ``25.0``.

    Assumptions
    -----------
    - Spherically symmetric layer geometry.
    - All values in MKS units.
    - Solidus/liquidus pressure dependence deferred to Phase 9.
    """

    def __cinit__(self, *args, **kwargs):
        self._solidliquid_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner_m,
            double radius_outer_m,
            double mass_kg,
            str    material_name                   = "",
            bint   is_tidal                        = True,
            double tidal_scale                     = 1.0,
            double shear_modulus_static_pa         = 0.0,
            double bulk_modulus_static_pa          = 0.0,
            double shear_viscosity_static_pas      = 0.0,
            double bulk_viscosity_static_pas       = 0.0,
            complex love_number_k                  = 0+0j,
            complex love_number_h                  = 0+0j,
            complex love_number_l                  = 0+0j,
            double thermal_conductivity_ref_w_mk   = 4.0,
            double thermal_expansion_ref_1_k       = 3.0e-5,
            double heat_capacity_ref_j_kgk         = 1200.0,
            double activation_energy_j_mol         = 300.0e3,
            double activation_volume_m3_mol        = 5.0e-6,
            double solidus_temperature_k           = 1600.0,
            double liquidus_temperature_k          = 2000.0,
            double melt_fraction_exponent          = 1.0,
            double reference_density_kg_m3         = 3500.0,
            double reference_temperature_k         = 1600.0,
            double melt_viscosity_reduction        = 25.0,
            str    tidal_scale_method              = "user_provided"):
        cdef c_SolidLiquidConfig config
        config.name                          = name.encode("utf-8")
        config.layer_index                   = layer_index
        config.radius_inner_m                = radius_inner_m
        config.radius_outer_m                = radius_outer_m
        config.mass_kg                       = mass_kg
        config.material_name                 = material_name.encode("utf-8")
        config.is_tidal                      = is_tidal
        config.tidal_scale                   = tidal_scale
        config.tidal_scale_method            = c_tidal_scale_method_from_name(tidal_scale_method.encode("utf-8"))
        config.shear_modulus_static_pa       = shear_modulus_static_pa
        config.bulk_modulus_static_pa        = bulk_modulus_static_pa
        config.shear_viscosity_static_pas    = shear_viscosity_static_pas
        config.bulk_viscosity_static_pas     = bulk_viscosity_static_pas
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.thermal_conductivity_ref_w_mk = thermal_conductivity_ref_w_mk
        config.thermal_expansion_ref_1_k     = thermal_expansion_ref_1_k
        config.heat_capacity_ref_j_kgk       = heat_capacity_ref_j_kgk
        config.activation_energy_j_mol       = activation_energy_j_mol
        config.activation_volume_m3_mol      = activation_volume_m3_mol
        config.solidus_temperature_k         = solidus_temperature_k
        config.liquidus_temperature_k        = liquidus_temperature_k
        config.melt_fraction_exponent        = melt_fraction_exponent
        config.reference_density_kg_m3       = reference_density_kg_m3
        config.reference_temperature_k       = reference_temperature_k
        config.melt_viscosity_reduction      = melt_viscosity_reduction
        cdef c_SolidLiquidLayer* raw = new c_SolidLiquidLayer(config)
        self._layer_ptr.reset(<c_BaseLayer*>raw)
        self._physics_ptr     = <c_PhysicsLayer*>raw
        self._solidliquid_ptr = raw
        self._ptr             = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._solidliquid_ptr = NULL  # base's unique_ptr owns the C++ object
        self._physics_ptr     = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Thermal property properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def thermal_conductivity_ref(self) -> float:
        """Reference thermal conductivity [W/(m·K)]."""
        return self._solidliquid_ptr.get_thermal_conductivity_ref()

    @property
    def thermal_expansion_ref(self) -> float:
        """Reference thermal expansion coefficient [1/K]."""
        return self._solidliquid_ptr.get_thermal_expansion_ref()

    @property
    def heat_capacity_ref(self) -> float:
        """Reference specific heat capacity [J/(kg·K)]."""
        return self._solidliquid_ptr.get_heat_capacity_ref()

    @property
    def activation_energy(self) -> float:
        """Arrhenius activation energy [J/mol]."""
        return self._solidliquid_ptr.get_activation_energy()

    @property
    def activation_volume(self) -> float:
        """Arrhenius activation volume [m³/mol]."""
        return self._solidliquid_ptr.get_activation_volume()

    @property
    def solidus_temperature(self) -> float:
        """Solidus temperature [K]."""
        return self._solidliquid_ptr.get_solidus_temperature()

    @property
    def liquidus_temperature(self) -> float:
        """Liquidus temperature [K]."""
        return self._solidliquid_ptr.get_liquidus_temperature()

    @property
    def melt_fraction_exponent(self) -> float:
        """Exponent in melt-fraction parameterization [dimensionless]."""
        return self._solidliquid_ptr.get_melt_fraction_exponent()

    @property
    def reference_density(self) -> float:
        """Reference density used for thermal diffusivity [kg/m³]."""
        return self._solidliquid_ptr.get_reference_density()

    @property
    def reference_temperature(self) -> float:
        """Reference temperature for Arrhenius viscosity [K]."""
        return self._solidliquid_ptr.get_reference_temperature()

    @property
    def melt_viscosity_reduction(self) -> float:
        """Exponential melt-viscosity reduction coefficient [dimensionless]."""
        return self._solidliquid_ptr.get_melt_viscosity_reduction()

    @property
    def cooling_set(self) -> bool:
        """True after a cooling sub-model has been attached."""
        return self._solidliquid_ptr.get_cooling_set()

    @property
    def radiogenics_set(self) -> bool:
        """True after a radiogenics sub-model has been attached."""
        return self._solidliquid_ptr.get_radiogenics_set()

    # ------------------------------------------------------------------------------------------------------------------
    # Sub-model attachment
    # ------------------------------------------------------------------------------------------------------------------
    def set_cooling(self, CoolingBase cooling not None):
        """Attach a cooling (heat-transport) sub-model.

        Ownership of the C++ model is transferred from ``cooling`` into this
        layer; the passed ``CoolingBase`` becomes an empty, non-owning shell and
        must not be reused.

        Parameters
        ----------
        cooling : CoolingBase
            A cooling model (e.g. ``make_cooling("convection")``).

        Raises
        ------
        ValueError
            If ``cooling`` has already been attached or otherwise moved.
        """
        if cooling._cooling_ptr.get() == NULL:
            raise ValueError(
                "This cooling model holds no C++ object (already attached or moved).")
        self._solidliquid_ptr.set_cooling(move(cooling._cooling_ptr))

    def set_radiogenics(self, RadiogenicsBase radiogenics not None):
        """Attach a radiogenic-heating sub-model.

        Ownership of the C++ model is transferred from ``radiogenics`` into this
        layer; the passed ``RadiogenicsBase`` becomes an empty, non-owning shell
        and must not be reused.

        Parameters
        ----------
        radiogenics : RadiogenicsBase
            A radiogenics model (e.g. ``make_radiogenics("fixed")``).

        Raises
        ------
        ValueError
            If ``radiogenics`` has already been attached or otherwise moved.
        """
        if radiogenics._radiogenics_ptr.get() == NULL:
            raise ValueError(
                "This radiogenics model holds no C++ object (already attached or moved).")
        self._solidliquid_ptr.set_radiogenics(move(radiogenics._radiogenics_ptr))

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_melt_fraction(self, double temperature_k, double pressure_pa = 0.0) -> float:
        """Volumetric melt fraction [0, 1] at temperature_k and pressure_pa.

        Uses a power-law interpolation between solidus and liquidus:
        φ = clamp((T - T_solidus)/(T_liquidus - T_solidus), 0, 1)^n

        Parameters
        ----------
        temperature_k : float
            Temperature [K].
        pressure_pa : float, optional
            Pressure [Pa]. Reserved for future pressure-dependent melt curve;
            currently unused. Default ``0.0``.

        Returns
        -------
        float
            Melt fraction [0, 1].
        """
        return self._solidliquid_ptr.calc_melt_fraction(temperature_k, pressure_pa)

    def calc_viscosity(self, double temperature_k, double pressure_pa = 0.0) -> float:
        """Effective viscosity [Pa·s] via Arrhenius + partial-melt reduction.

        η(T,P) = η_ref · exp((E_a + P·V_a)/(R·T) - E_a/(R·T_ref)) · exp(-C·φ)

        Parameters
        ----------
        temperature_k : float
            Temperature [K].
        pressure_pa : float, optional
            Pressure [Pa]. Default ``0.0``.

        Returns
        -------
        float
            Effective viscosity [Pa·s].
        """
        return self._solidliquid_ptr.calc_viscosity(temperature_k, pressure_pa)

    def calc_shear_modulus(self, double temperature_k, double pressure_pa = 0.0) -> float:
        """Effective shear modulus [Pa] accounting for melt fraction.

        G_eff = G_static · (1 - φ)

        Parameters
        ----------
        temperature_k : float
            Temperature [K].
        pressure_pa : float, optional
            Pressure [Pa]. Default ``0.0``.

        Returns
        -------
        float
            Effective shear modulus [Pa].
        """
        return self._solidliquid_ptr.calc_shear_modulus(temperature_k, pressure_pa)

    def calc_thermal_conductivity(self, double temperature_k) -> float:
        """Thermal conductivity [W/(m·K)].

        Returns the reference value (temperature dependence deferred to later phases).

        Parameters
        ----------
        temperature_k : float
            Temperature [K].

        Returns
        -------
        float
            Thermal conductivity [W/(m·K)].
        """
        return self._solidliquid_ptr.calc_thermal_conductivity(temperature_k)

    def calc_thermal_diffusivity(self, double temperature_k) -> float:
        """Thermal diffusivity [m²/s] = k / (ρ_ref · c_p).

        Parameters
        ----------
        temperature_k : float
            Temperature [K].

        Returns
        -------
        float
            Thermal diffusivity [m²/s].
        """
        return self._solidliquid_ptr.calc_thermal_diffusivity(temperature_k)

    def calc_adiabatic_temperature_gradient(self, double temperature_k,
                                            double pressure_pa = 0.0) -> float:
        """Adiabatic temperature gradient [K/m] = α · T · g / c_p.

        Gravity is taken from the EOS profile at the outer boundary if
        available; returns 0.0 when EOS data has not been populated.

        Parameters
        ----------
        temperature_k : float
            Temperature [K].
        pressure_pa : float, optional
            Pressure [Pa] (reserved; currently unused). Default ``0.0``.

        Returns
        -------
        float
            Adiabatic temperature gradient [K/m].
        """
        return self._solidliquid_ptr.calc_adiabatic_temperature_gradient(
            temperature_k, pressure_pa)

    def calc_heat_flux_conductive(self, double temperature_base_k,
                                  double temperature_top_k) -> float:
        """Conductive heat flux [W/m²] = k · (T_base - T_top) / thickness.

        Parameters
        ----------
        temperature_base_k : float
            Temperature at the base (inner boundary) [K].
        temperature_top_k : float
            Temperature at the top (outer boundary) [K].

        Returns
        -------
        float
            Conductive heat flux [W/m²].  Positive when T_base > T_top.
        """
        return self._solidliquid_ptr.calc_heat_flux_conductive(
            temperature_base_k, temperature_top_k)

    def calc_radiogenic_heating(self, double time_s, double mass_kg) -> float:
        """Radiogenic heating [W] from the attached sub-model.

        Returns 0.0 when no radiogenics sub-model has been attached (Phase 7).

        Parameters
        ----------
        time_s : float
            Elapsed time since reference epoch [s].
        mass_kg : float
            Radiogenic mass [kg].

        Returns
        -------
        float
            Radiogenic heating power [W].
        """
        return self._solidliquid_ptr.calc_radiogenic_heating(time_s, mass_kg)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            All BaseLayer + PhysicsLayer keys plus the 11 SolidLiquidLayer
            thermal/melt parameters.
        """
        d = PhysicsLayer.get_config_dict(self)
        d["thermal_conductivity_ref_w_mk"] = self._solidliquid_ptr.get_thermal_conductivity_ref()
        d["thermal_expansion_ref_1_k"]     = self._solidliquid_ptr.get_thermal_expansion_ref()
        d["heat_capacity_ref_j_kgk"]       = self._solidliquid_ptr.get_heat_capacity_ref()
        d["activation_energy_j_mol"]       = self._solidliquid_ptr.get_activation_energy()
        d["activation_volume_m3_mol"]      = self._solidliquid_ptr.get_activation_volume()
        d["solidus_temperature_k"]         = self._solidliquid_ptr.get_solidus_temperature()
        d["liquidus_temperature_k"]        = self._solidliquid_ptr.get_liquidus_temperature()
        d["melt_fraction_exponent"]        = self._solidliquid_ptr.get_melt_fraction_exponent()
        d["reference_density_kg_m3"]       = self._solidliquid_ptr.get_reference_density()
        d["reference_temperature_k"]       = self._solidliquid_ptr.get_reference_temperature()
        d["melt_viscosity_reduction"]      = self._solidliquid_ptr.get_melt_viscosity_reduction()
        return d
