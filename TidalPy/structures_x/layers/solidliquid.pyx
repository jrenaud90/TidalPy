# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's solid/liquid layer class.

SolidLiquidLayer extends PhysicsLayer with optional cooling and radiogenics sub-models and the conductive and
adiabatic calculations that need the layer's geometry or solved profile.
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_unique

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass, c_PhysicsBase, cy_physics_model_config
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer, c_tidal_scale_method_from_name
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers
from TidalPy.cooling_x.cooling cimport CoolingBase
from TidalPy.radiogenics_x.radiogenics cimport RadiogenicsBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class SolidLiquidLayer(PhysicsLayer):
    """Thermo-mechanical layer with optional cooling and radiogenics sub-models.

    Extends PhysicsLayer with thermal diffusivity, the adiabatic gradient, and conductive heat flux, which read
    the thermal constants (conductivity, expansivity, heat capacity) of the layer's material, its EOS model.
    Viscosity, melt fraction, and the melt-reduced shear modulus belong to the material too; the EOS solve
    evaluates them, and the radius getters read them back.

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
    is_volume_fixed : bool, optional
        False lets the layer grow or shrink to hold its mass during an EOS solve. Default ``True``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        Dimensionless tidal heating scale. Default ``1.0``.
    love_number_k : complex, optional
        Potential Love number k (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number h (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number l (placeholder). Default ``0+0j``.
    tidal_scale_method : str, optional
        How the layer's share of the world's tidal heating is set. Default ``"user_provided"``.
    is_solid : bool, optional
        False marks the layer liquid for the radial Love-number solver. Default ``True``.
    is_static : bool, optional
        Use the static (no inertia) approximation in the radial solver. Default ``True``, so a liquid layer
        is a static liquid unless this is set False.
    is_incompressible : bool, optional
        Use the incompressible approximation in the radial solver. Default ``False``.
    temperature : float, optional
        Layer temperature [K] at which its viscosity and melt models are evaluated. Default ``0.0``, the cold
        rigid limit of the viscosity laws.
    use_thermal_eos : bool, optional
        Pass the temperature to the EOS model, so the density and bulk modulus depend on it. Default ``False``.
    use_heating : bool, optional
        Let the world's heat sources (this layer's radiogenics model among them) act inside the layer during a
        thermal EOS solve. Default ``False``.

    Assumptions
    -----------
    - Spherically symmetric layer geometry.
    """

    def __cinit__(self, *args, **kwargs):
        self._solidliquid_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner,
            double radius_outer,
            double mass,
            str    material_name            = "",
            cpp_bool is_tidal               = True,
            cpp_bool is_volume_fixed        = True,
            double tidal_scale              = 1.0,
            complex love_number_k           = 0+0j,
            complex love_number_h           = 0+0j,
            complex love_number_l           = 0+0j,
            str    tidal_scale_method       = "user_provided",
            cpp_bool   is_solid             = True,
            cpp_bool   is_static            = True,
            cpp_bool   is_incompressible    = False,
            double temperature           = 0.0,
            cpp_bool use_thermal_eos     = False,
            cpp_bool use_heating         = False):
        cdef c_SolidLiquidConfig config
        config.name                 = name.encode("utf-8")
        config.layer_index          = layer_index
        config.radius_inner         = radius_inner
        config.radius_outer         = radius_outer
        config.mass                 = mass
        config.material_name        = material_name.encode("utf-8")
        config.is_tidal             = is_tidal
        config.is_volume_fixed      = is_volume_fixed
        config.tidal_scale          = tidal_scale
        config.tidal_scale_method   = c_tidal_scale_method_from_name(tidal_scale_method.encode("utf-8"))
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.is_solid             = is_solid
        config.is_static            = is_static
        config.is_incompressible    = is_incompressible
        config.temperature       = temperature
        config.use_thermal_eos   = use_thermal_eos
        config.use_heating       = use_heating
        # make_unique owns the allocation; ownership then moves into the base-typed member
        # (Cython cannot assign a unique_ptr[Derived] to a unique_ptr[Base] directly).
        cdef unique_ptr[c_SolidLiquidLayer] built = make_unique[c_SolidLiquidLayer](config)
        self._solidliquid_ptr = built.get()
        self._physics_ptr     = <c_PhysicsLayer*>self._solidliquid_ptr
        self._layer_ptr.reset(<c_BaseLayer*>built.release())
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()

    def __dealloc__(self):
        self._solidliquid_ptr = NULL  # base's unique_ptr owns the C++ object
        self._physics_ptr     = NULL

    @staticmethod
    cdef SolidLiquidLayer _view(c_SolidLiquidLayer* ptr, object world):
        cdef SolidLiquidLayer v = SolidLiquidLayer.__new__(SolidLiquidLayer)
        v._solidliquid_ptr = ptr
        v._physics_ptr     = <c_PhysicsLayer*>ptr
        v._init_view(<c_BaseLayer*>ptr, world)
        return v

    @property
    def thermal_conductivity(self) -> float:
        """Thermal conductivity k [W/(m K)] of the layer's material (its EOS model); NaN when none is attached."""
        return self._solidliquid_ptr.get_thermal_conductivity()

    @property
    def thermal_expansion(self) -> float:
        """Thermal expansivity alpha [1/K] of the layer's material; NaN when none is attached."""
        return self._solidliquid_ptr.get_thermal_expansion()

    @property
    def heat_capacity(self) -> float:
        """Specific heat capacity c_p [J/(kg K)] of the layer's material; NaN when none is attached."""
        return self._solidliquid_ptr.get_heat_capacity()

    @property
    def cooling_set(self) -> bool:
        """True after a cooling sub-model has been attached."""
        return self._solidliquid_ptr.get_cooling_set()

    @property
    def radiogenics_set(self) -> bool:
        """True after a radiogenics sub-model has been attached."""
        return self._solidliquid_ptr.get_radiogenics_set()

    def set_cooling(self, CoolingBase cooling not None):
        """Attach a cooling (heat-transport) sub-model.

        Ownership of the C++ model moves out of ``cooling``, which is left an empty shell and must not be reused.

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
        cooling._ptr = NULL

    def set_radiogenics(self, RadiogenicsBase radiogenics not None):
        """Attach a radiogenic-heating sub-model.

        Ownership of the C++ model moves out of ``radiogenics``, which is left an empty shell and must not be
        reused.

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
        radiogenics._ptr = NULL

    def calc_thermal_conductivity(self, double temperature) -> float:
        """Thermal conductivity [W/(m·K)]: the reference value, with no temperature dependence modeled.

        Parameters
        ----------
        temperature : float
            Temperature [K].

        Returns
        -------
        float
            Thermal conductivity [W/(m·K)].
        """
        return self._solidliquid_ptr.calc_thermal_conductivity(temperature)

    def calc_thermal_diffusivity(self, double temperature) -> float:
        """Thermal diffusivity [m²/s] = k / (ρ_ref · c_p).

        Parameters
        ----------
        temperature : float
            Temperature [K].

        Returns
        -------
        float
            Thermal diffusivity [m²/s].
        """
        return self._solidliquid_ptr.calc_thermal_diffusivity(temperature)

    def calc_adiabatic_temperature_gradient(self, double temperature,
                                            double pressure = 0.0) -> float:
        """Adiabatic temperature gradient [K/m] = α · T · g / c_p.

        Gravity comes from the EOS profile at the outer boundary; 0.0 when that profile is unpopulated.

        Parameters
        ----------
        temperature : float
            Temperature [K].
        pressure : float, optional
            Pressure [Pa]. Not used. Default ``0.0``.

        Returns
        -------
        float
            Adiabatic temperature gradient [K/m].
        """
        return self._solidliquid_ptr.calc_adiabatic_temperature_gradient(
            temperature, pressure)

    def calc_heat_flux_conductive(self, double temperature_base,
                                  double temperature_top) -> float:
        """Conductive heat flux [W/m²] = k · (T_base - T_top) / thickness.

        Parameters
        ----------
        temperature_base : float
            Temperature at the base (inner boundary) [K].
        temperature_top : float
            Temperature at the top (outer boundary) [K].

        Returns
        -------
        float
            Conductive heat flux [W/m²]. Positive when T_base > T_top.
        """
        return self._solidliquid_ptr.calc_heat_flux_conductive(
            temperature_base, temperature_top)

    def calc_radiogenic_heating(self, double time, double mass) -> float:
        """Radiogenic heating [W] from the attached sub-model; 0.0 when none is attached.

        Parameters
        ----------
        time : float
            Elapsed time since reference epoch [s].
        mass : float
            Radiogenic mass [kg].

        Returns
        -------
        float
            Radiogenic heating power [W].
        """
        return self._solidliquid_ptr.calc_radiogenic_heating(time, mass)

    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            The PhysicsLayer keys plus the ``cooling`` and ``radiogenics`` sub-tables when those models are
            attached. The thermal constants are the material's, so they sit in the ``material`` table.
        """
        cdef dict d = PhysicsLayer.get_config_dict(self)
        cdef const c_PhysicsBase* model_ptr
        model_ptr = <const c_PhysicsBase*>self._solidliquid_ptr.get_cooling_model()
        if model_ptr != NULL:
            d["cooling"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>self._solidliquid_ptr.get_radiogenics_model()
        if model_ptr != NULL:
            d["radiogenics"] = cy_physics_model_config(model_ptr)
        return d
