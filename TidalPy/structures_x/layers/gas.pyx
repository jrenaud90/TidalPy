# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
gas.pyx
Cython/Python wrapper for TidalPy's gas layer class.

GasLayer: extends PhysicsLayer with ideal-gas thermodynamic properties
(adiabatic lapse rate, scale height, ideal-gas pressure, sound speed).
No phase changes, no solidus/liquidus, and no cooling or radiogenics
sub-models.
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# GasLayer
# =====================================================================================================================
cdef class GasLayer(PhysicsLayer):
    """Ideal-gas/fluid layer with thermodynamic properties.

    Extends PhysicsLayer with ideal-gas calculations:

    - Dry adiabatic lapse rate.
    - Barometric scale height.
    - Ideal-gas pressure from density and temperature.
    - Adiabatic sound speed.

    No phase changes, cooling, or radiogenics sub-models are available
    (use SolidLiquidLayer for those features).

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
        Unrelaxed shear modulus [Pa]. Default ``0.0``.
    bulk_modulus_static_pa : float, optional
        Unrelaxed bulk modulus [Pa]. Default ``0.0``.
    shear_viscosity_static_pas : float, optional
        Reference dynamic shear viscosity [Pa·s]. Default ``0.0``.
    bulk_viscosity_static_pas : float, optional
        Reference dynamic bulk viscosity [Pa·s]. Default ``0.0``.
    love_number_k : complex, optional
        Potential Love number k (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number h (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number l (placeholder). Default ``0+0j``.
    mean_molecular_weight_kg_mol : float, optional
        Mean molecular weight of the gas [kg/mol]. Default ``2e-3`` (H₂).
    adiabatic_index : float, optional
        Ratio of specific heats γ = c_p/c_v [dimensionless]. Default ``1.4``.
    reference_temperature_k : float, optional
        Reference temperature [K]. Default ``300.0``.
    reference_density_kg_m3 : float, optional
        Reference density [kg/m³]. Default ``1.0``.

    Assumptions
    -----------
    - Spherically symmetric layer geometry.
    - All values in MKS units.
    - Gas is ideal; real-gas corrections are not included.
    - Universal gas constant sourced from TidalPy global configuration.
    """

    def __cinit__(self, *args, **kwargs):
        self._gas_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner_m,
            double radius_outer_m,
            double mass_kg,
            str    material_name                = "",
            bint   is_tidal                     = True,
            double tidal_scale                  = 1.0,
            double shear_modulus_static_pa      = 0.0,
            double bulk_modulus_static_pa       = 0.0,
            double shear_viscosity_static_pas   = 0.0,
            double bulk_viscosity_static_pas    = 0.0,
            complex love_number_k               = 0+0j,
            complex love_number_h               = 0+0j,
            complex love_number_l               = 0+0j,
            double mean_molecular_weight_kg_mol = 2.0e-3,
            double adiabatic_index              = 1.4,
            double reference_temperature_k      = 300.0,
            double reference_density_kg_m3      = 1.0):
        cdef c_GasConfig config
        config.name                          = name.encode("utf-8")
        config.layer_index                   = layer_index
        config.radius_inner_m                = radius_inner_m
        config.radius_outer_m                = radius_outer_m
        config.mass_kg                       = mass_kg
        config.material_name                 = material_name.encode("utf-8")
        config.is_tidal                      = is_tidal
        config.tidal_scale                   = tidal_scale
        config.shear_modulus_static_pa       = shear_modulus_static_pa
        config.bulk_modulus_static_pa        = bulk_modulus_static_pa
        config.shear_viscosity_static_pas    = shear_viscosity_static_pas
        config.bulk_viscosity_static_pas     = bulk_viscosity_static_pas
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.mean_molecular_weight_kg_mol  = mean_molecular_weight_kg_mol
        config.adiabatic_index               = adiabatic_index
        config.reference_temperature_k       = reference_temperature_k
        config.reference_density_kg_m3       = reference_density_kg_m3
        cdef c_GasLayer* raw = new c_GasLayer(config)
        self._layer_ptr.reset(<c_BaseLayer*>raw)
        self._physics_ptr = <c_PhysicsLayer*>raw
        self._gas_ptr     = raw
        self._ptr         = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._gas_ptr     = NULL  # base's unique_ptr owns the C++ object
        self._physics_ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Gas property properties
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
    def calc_adiabatic_lapse_rate(self, double gravity_m_s2) -> float:
        """Dry adiabatic lapse rate [K/m] = g * (γ-1) * M / (γ * R).

        Parameters
        ----------
        gravity_m_s2 : float
            Gravitational acceleration [m/s²].

        Returns
        -------
        float
            Adiabatic lapse rate [K/m].  Returns 0.0 for invalid inputs.

        Assumptions
        -----------
        - Ideal gas.
        - R from TidalPy global configuration.
        """
        return self._gas_ptr.calc_adiabatic_lapse_rate(gravity_m_s2)

    def calc_scale_height(self, double temperature_k,
                          double gravity_m_s2) -> float:
        """Barometric scale height [m] = R * T / (g * M).

        Parameters
        ----------
        temperature_k : float
            Temperature [K].
        gravity_m_s2 : float
            Gravitational acceleration [m/s²].

        Returns
        -------
        float
            Scale height [m].  Returns 0.0 for invalid inputs.

        Assumptions
        -----------
        - Ideal gas.
        - R from TidalPy global configuration.
        """
        return self._gas_ptr.calc_scale_height(temperature_k, gravity_m_s2)

    def calc_pressure_ideal_gas(self, double temperature_k,
                                double density_kg_m3) -> float:
        """Ideal-gas pressure [Pa] = ρ * R * T / M.

        Parameters
        ----------
        temperature_k : float
            Temperature [K].
        density_kg_m3 : float
            Gas density [kg/m³].

        Returns
        -------
        float
            Pressure [Pa].  Returns 0.0 for invalid inputs.

        Assumptions
        -----------
        - Ideal gas.
        - R from TidalPy global configuration.
        """
        return self._gas_ptr.calc_pressure_ideal_gas(temperature_k, density_kg_m3)

    def calc_sound_speed(self, double temperature_k) -> float:
        """Adiabatic sound speed [m/s] = sqrt(γ * R * T / M).

        Parameters
        ----------
        temperature_k : float
            Temperature [K].

        Returns
        -------
        float
            Sound speed [m/s].  Returns 0.0 for invalid inputs.

        Assumptions
        -----------
        - Ideal gas.
        - R from TidalPy global configuration.
        """
        return self._gas_ptr.calc_sound_speed(temperature_k)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            All BaseLayer + PhysicsLayer keys plus the 4 GasLayer
            thermodynamic parameters.
        """
        d = PhysicsLayer.get_config_dict(self)
        d["mean_molecular_weight_kg_mol"] = self._gas_ptr.get_mean_molecular_weight()
        d["adiabatic_index"]              = self._gas_ptr.get_adiabatic_index()
        d["reference_temperature_k"]      = self._gas_ptr.get_reference_temperature()
        d["reference_density_kg_m3"]      = self._gas_ptr.get_reference_density()
        return d
