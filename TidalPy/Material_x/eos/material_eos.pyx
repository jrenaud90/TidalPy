# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for the material EOS models.

A model is a layer's material. It returns the density [kg/m^3] from the local pressure [Pa] (analytic models) or
radius [m] (interpolated model), and it owns every other frequency-independent property: the static shear law,
the static bulk modulus and viscosities, and the optional viscosity and partial-melt models. The whole-planet EOS
solve evaluates all of it as it integrates. Models: ConstantDensityEOS ("constant"), BirchMurnaghanEOS ("bm"),
VinetEOS ("vinet"), InterpolatedEOS ("interp").

Every model takes a thermal expansivity and a reference temperature. Birch-Murnaghan and Vinet add the thermal
pressure alpha0 K0 (T - T_ref) to their pressure law; the constant and interpolated models scale their density by
exp(-alpha0 (T - T_ref)). A zero expansivity (the default) or no temperature gives the athermal EOS.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport (
    PhysicsBase, c_TidalPyBaseClass, c_PhysicsBase, cy_physics_model_config)
from TidalPy.Utilities_x.classes_x.classes import check_config_keys, factory_defaults
from TidalPy.viscosity_x.viscosity cimport ViscosityBase
from TidalPy.partial_melt_x.partial_melt cimport PartialMeltBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# Free analytic pressure laws (Python-accessible for cross-checks)
# =====================================================================================================================
def birch_murnaghan_pressure(double eta, double reference_bulk_modulus,
                             double bulk_modulus_derivative=4.0) -> float:
    """3rd-order Birch-Murnaghan pressure [Pa] at compression eta = rho/rho0."""
    return eos_bm_pressure(eta, reference_bulk_modulus, bulk_modulus_derivative)


def vinet_pressure(double eta, double reference_bulk_modulus,
                   double bulk_modulus_derivative=4.0) -> float:
    """Vinet pressure [Pa] at compression eta = rho/rho0."""
    return eos_vinet_pressure(eta, reference_bulk_modulus, bulk_modulus_derivative)


# =====================================================================================================================
# MaterialEOSBase
# =====================================================================================================================
cdef class MaterialEOSBase(PhysicsBase):
    """Abstract base for material EOS models. Instantiate a concrete subclass."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_MaterialEOSBase> auto-inits to nullptr

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "MaterialEOSBase is abstract; instantiate a concrete model "
            "(ConstantDensityEOS, BirchMurnaghanEOS, VinetEOS, InterpolatedEOS).")

    def __dealloc__(self):
        self._eos_ptr.reset()
        self._ptr = NULL

    cdef c_MaterialEOSBase* _model(self) except NULL:
        # Attaching the model to a layer moves it out of this wrapper, which is then an empty shell.
        cdef c_MaterialEOSBase* model = self._eos_ptr.get()
        if model == NULL:
            raise ValueError("This EOS model holds no C++ object (already attached or moved).")
        return model

    def calc_density(self, double pressure, temperature=None, double radius=0.0) -> float:
        """Density [kg/m^3] at pressure [Pa] (analytic models) or radius [m] (interpolated model).

        ``temperature`` [K] enters through the model's thermal expansivity; ``None`` gives the athermal density.
        """
        cdef double temperature_value = d_NAN if temperature is None else <double>temperature
        return self._model().calc_density(pressure, temperature_value, radius)

    def calc_bulk_modulus(self, double pressure, temperature=None, double radius=0.0) -> float:
        """Isothermal bulk modulus K = rho dP/drho [Pa] at pressure [Pa] and temperature [K].

        Birch-Murnaghan and Vinet return the modulus of their pressure law at the solved compression. The other
        models return their bulk table at ``radius`` [m], or NaN when they have none.
        """
        cdef double temperature_value = d_NAN if temperature is None else <double>temperature
        return self._model().calc_bulk_modulus(pressure, temperature_value, radius)

    @property
    def thermal_expansion(self) -> float:
        """Thermal expansivity alpha0 [1/K]; zero is the athermal EOS."""
        return self._model().get_thermal_expansion()

    @property
    def reference_temperature(self) -> float:
        """Temperature [K] at which the reference density (and bulk modulus) apply."""
        return self._model().get_reference_temperature()

    def get_tabulated_shear_modulus(self, double radius) -> float:
        """Static shear modulus [Pa] at radius [m] from the model's table; NaN unless it carries one."""
        return self._model().get_tabulated_shear_modulus(radius)

    def get_tabulated_bulk_modulus(self, double radius) -> float:
        """Static bulk modulus [Pa] at radius [m] from the model's table; NaN unless it carries one."""
        return self._model().get_tabulated_bulk_modulus(radius)

    def get_tabulated_shear_viscosity(self, double radius) -> float:
        """Shear viscosity [Pa s] at radius [m] from the model's table; NaN unless it carries one."""
        return self._model().get_tabulated_shear_viscosity(radius)

    def get_tabulated_bulk_viscosity(self, double radius) -> float:
        """Bulk viscosity [Pa s] at radius [m] from the model's table; NaN unless it carries one."""
        return self._model().get_tabulated_bulk_viscosity(radius)

    # ------------------------------------------------------------------------------------------------------------------
    # The material: static constants, the shear law, and the viscosity and partial-melt models
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def shear_modulus_static(self) -> float:
        """Static (unrelaxed) shear modulus mu0 [Pa] of the shear law."""
        return self._model().get_shear_modulus_static()

    @shear_modulus_static.setter
    def shear_modulus_static(self, double value):
        self._model().set_shear_modulus_static(value)

    @property
    def bulk_modulus_static(self) -> float:
        """Static bulk modulus [Pa], used by a model with no pressure law or bulk table of its own."""
        return self._model().get_bulk_modulus_static()

    @bulk_modulus_static.setter
    def bulk_modulus_static(self, double value):
        self._model().set_bulk_modulus_static(value)

    @property
    def shear_viscosity_static(self) -> float:
        """Static shear viscosity [Pa s], used when no shear viscosity model is attached; NaN until set."""
        return self._model().get_shear_viscosity_static()

    @shear_viscosity_static.setter
    def shear_viscosity_static(self, double value):
        self._model().set_shear_viscosity_static(value)

    @property
    def bulk_viscosity_static(self) -> float:
        """Static bulk viscosity [Pa s], used when no bulk viscosity model is attached; NaN until set."""
        return self._model().get_bulk_viscosity_static()

    @bulk_viscosity_static.setter
    def bulk_viscosity_static(self, double value):
        self._model().set_bulk_viscosity_static(value)

    @property
    def shear_modulus_pressure_derivative(self) -> float:
        """Pressure derivative of the static shear law [dimensionless]."""
        return self._model().get_shear_modulus_pressure_derivative()

    @property
    def shear_modulus_temperature_derivative(self) -> float:
        """Temperature derivative of the static shear law [Pa/K]."""
        return self._model().get_shear_modulus_temperature_derivative()

    @property
    def shear_modulus_reference_temperature(self) -> float:
        """Reference temperature of the static shear law [K]."""
        return self._model().get_shear_modulus_reference_temperature()

    @property
    def thermal_conductivity(self) -> float:
        """Thermal conductivity k [W/(m K)]."""
        return self._model().get_thermal_conductivity()

    @property
    def heat_capacity(self) -> float:
        """Specific heat capacity c_p [J/(kg K)]."""
        return self._model().get_heat_capacity()

    def calc_thermal_diffusivity(self, double density) -> float:
        """Thermal diffusivity k / (rho c_p) [m^2/s] at a density [kg/m^3]."""
        return self._model().calc_thermal_diffusivity(density)

    def set_shear_viscosity(self, ViscosityBase viscosity not None):
        """Attach a viscosity model supplying the shear viscosity before the partial-melt model.

        Ownership of the C++ model moves out of ``viscosity``, which is left an empty shell and must not be reused.
        """
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError("This viscosity model holds no C++ object (already attached or moved).")
        self._model().set_shear_viscosity(move(viscosity._visc_ptr))

    def set_bulk_viscosity(self, ViscosityBase viscosity not None):
        """Attach a viscosity model supplying the bulk viscosity before the partial-melt model.

        Ownership of the C++ model moves out of ``viscosity``, which is left an empty shell and must not be reused.
        """
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError("This viscosity model holds no C++ object (already attached or moved).")
        self._model().set_bulk_viscosity(move(viscosity._visc_ptr))

    def set_partial_melt(self, PartialMeltBase partial_melt not None):
        """Attach a partial-melt model that weakens the static moduli and viscosities with melt fraction.

        Ownership of the C++ model moves out of ``partial_melt``, which is left an empty shell and must not be
        reused.
        """
        if partial_melt._melt_ptr.get() == NULL:
            raise ValueError("This partial-melt model holds no C++ object (already attached or moved).")
        self._model().set_partial_melt(move(partial_melt._melt_ptr))

    @property
    def shear_viscosity_set(self) -> bool:
        """True once a shear viscosity model is attached."""
        return self._model().get_shear_viscosity_model() != NULL

    @property
    def bulk_viscosity_set(self) -> bool:
        """True once a bulk viscosity model is attached."""
        return self._model().get_bulk_viscosity_model() != NULL

    @property
    def partial_melt_set(self) -> bool:
        """True once a partial-melt model is attached."""
        return self._model().get_partial_melt_model() != NULL

    def calc_material_state(self, double pressure, temperature=None, double radius=0.0,
                            thermal_density=True) -> dict:
        """Every frequency-independent property of the material at one point.

        This is what the EOS solve evaluates as it integrates: the shear law and constants, then a viscosity
        model in place of a constant, then a table or the pressure law's own bulk modulus in place of either,
        then the partial-melt model. After a solve, read these from the world or layer getters instead, which
        report what the solve used.

        Parameters
        ----------
        pressure : float
            Pressure [Pa].
        temperature : float, optional
            Temperature [K]. ``None`` is the cold limit of the viscosity laws and the athermal density.
        radius : float, optional
            Radius [m], read by the interpolated model alone.
        thermal_density : bool, optional
            Whether the density law sees the temperature (the viscosity and melt models always do).

        Returns
        -------
        dict
            ``density`` [kg/m^3], ``melt_fraction``, ``shear_modulus`` and ``bulk_modulus`` [Pa],
            ``shear_viscosity`` and ``bulk_viscosity`` [Pa s], all after the partial-melt model.
        """
        cdef double temperature_value = d_NAN if temperature is None else <double>temperature
        cdef c_MaterialState state
        self._model().calc_material_state(pressure, temperature_value, bool(thermal_density), radius, state)
        return {
            "density":         state.density,
            "melt_fraction":   state.melt_fraction,
            "shear_modulus":   state.shear_modulus,
            "bulk_modulus":    state.bulk_modulus,
            "shear_viscosity": state.shear_viscosity,
            "bulk_viscosity":  state.bulk_viscosity,
        }

    cpdef dict get_config_dict(self):
        """The model's configuration, with a sub-table for each attached viscosity or partial-melt model."""
        return cy_material_config(self._model())


cdef dict cy_material_config(const c_MaterialEOSBase* eos_ptr):
    cdef dict d = cy_physics_model_config(<const c_PhysicsBase*>eos_ptr)
    cdef const c_PhysicsBase* model_ptr
    model_ptr = <const c_PhysicsBase*>eos_ptr.get_shear_viscosity_model()
    if model_ptr != NULL:
        d["shear_viscosity"] = cy_physics_model_config(model_ptr)
    model_ptr = <const c_PhysicsBase*>eos_ptr.get_bulk_viscosity_model()
    if model_ptr != NULL:
        d["bulk_viscosity"] = cy_physics_model_config(model_ptr)
    model_ptr = <const c_PhysicsBase*>eos_ptr.get_partial_melt_model()
    if model_ptr != NULL:
        d["partial_melt"] = cy_physics_model_config(model_ptr)
    return d


# The keyword arguments every model constructor takes for the material, mapped onto c_MaterialEOSConfig.
_MATERIAL_KWARGS = (
    "shear_modulus_static", "bulk_modulus_static", "shear_viscosity_static", "bulk_viscosity_static",
    "shear_modulus_pressure_derivative", "shear_modulus_temperature_derivative",
    "shear_modulus_reference_temperature", "thermal_conductivity", "heat_capacity")


cdef int cy_fill_material_config(c_MaterialEOSConfig& config, dict material) except -1:
    for key in material:
        if key not in _MATERIAL_KWARGS:
            raise TypeError(f"unexpected material keyword argument {key!r}; expected one of {_MATERIAL_KWARGS}")
    if "shear_modulus_static" in material:
        config.shear_modulus_static = material["shear_modulus_static"]
    if "bulk_modulus_static" in material:
        config.bulk_modulus_static = material["bulk_modulus_static"]
    if "shear_viscosity_static" in material:
        config.shear_viscosity_static = material["shear_viscosity_static"]
    if "bulk_viscosity_static" in material:
        config.bulk_viscosity_static = material["bulk_viscosity_static"]
    if "shear_modulus_pressure_derivative" in material:
        config.shear_modulus_pressure_derivative = material["shear_modulus_pressure_derivative"]
    if "shear_modulus_temperature_derivative" in material:
        config.shear_modulus_temperature_derivative = material["shear_modulus_temperature_derivative"]
    if material.get("shear_modulus_reference_temperature") is not None:
        config.shear_modulus_reference_temperature = material["shear_modulus_reference_temperature"]
    if "thermal_conductivity" in material:
        config.thermal_conductivity = material["thermal_conductivity"]
    if "heat_capacity" in material:
        config.heat_capacity = material["heat_capacity"]
    return 0


# =====================================================================================================================
# EOS Models
# =====================================================================================================================
cdef class ConstantDensityEOS(MaterialEOSBase):
    """Incompressible (uniform) density EOS."""

    def __cinit__(self, *args, **kwargs):
        self._constant_ptr = NULL

    def __init__(self, double reference_density=3500.0, double thermal_expansion=0.0, reference_temperature=None,
                 **material):
        # None keeps the C++ default reference temperature. The material keywords are those of
        # _MATERIAL_KWARGS (static moduli and viscosities, and the shear law).
        cdef c_MaterialEOSConfig config
        cy_fill_material_config(config, material)
        config.reference_density = reference_density
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.Constant, config)
        self._constant_ptr = <c_ConstantDensityEOS*>ptr.get()
        self._eos_ptr      = move(ptr)
        self._ptr          = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._constant_ptr = NULL

    @property
    def reference_density(self) -> float:
        """Reference (uniform) density [kg/m^3]."""
        return self._constant_ptr.get_reference_density()


cdef class BirchMurnaghanEOS(MaterialEOSBase):
    """3rd-order Birch-Murnaghan EOS; density from pressure."""

    def __cinit__(self, *args, **kwargs):
        self._bm_ptr = NULL

    def __init__(
            self,
            double reference_density=3500.0,
            double reference_bulk_modulus=1.0e11,
            double bulk_modulus_derivative=4.0,
            invert_rtol=None,
            invert_max_iters=None,
            double thermal_expansion=0.0,
            reference_temperature=None,
            **material):
        # None keeps the C++ default inversion settings and reference temperature.
        cdef c_MaterialEOSConfig config
        cy_fill_material_config(config, material)
        config.reference_density   = reference_density
        config.reference_bulk_modulus = reference_bulk_modulus
        config.bulk_modulus_derivative   = bulk_modulus_derivative
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        if invert_rtol is not None:
            config.invert_rtol = <double>invert_rtol
        if invert_max_iters is not None:
            config.invert_max_iters = <int>invert_max_iters
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.BirchMurnaghan, config)
        self._bm_ptr  = <c_BirchMurnaghanEOS*>ptr.get()
        self._eos_ptr = move(ptr)
        self._ptr     = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._bm_ptr = NULL

    @property
    def reference_density(self) -> float:
        """Reference density rho0 [kg/m^3]."""
        return self._bm_ptr.get_reference_density()

    @property
    def reference_bulk_modulus(self) -> float:
        """Reference bulk modulus K0 [Pa]."""
        return self._bm_ptr.get_reference_bulk_modulus()

    @property
    def bulk_modulus_derivative(self) -> float:
        """Pressure derivative of the bulk modulus K0' [dimensionless]."""
        return self._bm_ptr.get_bulk_modulus_derivative()

    @property
    def invert_rtol(self) -> float:
        """Relative convergence tolerance for the density-from-pressure inversion."""
        return self._bm_ptr.get_invert_rtol()

    @property
    def invert_max_iters(self) -> int:
        """Hard iteration cap (termination safeguard) for the inversion."""
        return self._bm_ptr.get_invert_max_iters()


cdef class VinetEOS(MaterialEOSBase):
    """Vinet (universal) EOS; density from pressure."""

    def __cinit__(self, *args, **kwargs):
        self._vinet_ptr = NULL

    def __init__(
            self,
            double reference_density=3500.0,
            double reference_bulk_modulus=1.0e11,
            double bulk_modulus_derivative=4.0,
            invert_rtol=None,
            invert_max_iters=None,
            double thermal_expansion=0.0,
            reference_temperature=None,
            **material):
        # None keeps the C++ default inversion settings and reference temperature.
        cdef c_MaterialEOSConfig config
        cy_fill_material_config(config, material)
        config.reference_density   = reference_density
        config.reference_bulk_modulus = reference_bulk_modulus
        config.bulk_modulus_derivative   = bulk_modulus_derivative
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        if invert_rtol is not None:
            config.invert_rtol = <double>invert_rtol
        if invert_max_iters is not None:
            config.invert_max_iters = <int>invert_max_iters
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.Vinet, config)
        self._vinet_ptr = <c_VinetEOS*>ptr.get()
        self._eos_ptr   = move(ptr)
        self._ptr       = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._vinet_ptr = NULL

    @property
    def reference_density(self) -> float:
        """Reference density rho0 [kg/m^3]."""
        return self._vinet_ptr.get_reference_density()

    @property
    def reference_bulk_modulus(self) -> float:
        """Reference bulk modulus K0 [Pa]."""
        return self._vinet_ptr.get_reference_bulk_modulus()

    @property
    def bulk_modulus_derivative(self) -> float:
        """Pressure derivative of the bulk modulus K0' [dimensionless]."""
        return self._vinet_ptr.get_bulk_modulus_derivative()

    @property
    def invert_rtol(self) -> float:
        """Relative convergence tolerance for the density-from-pressure inversion."""
        return self._vinet_ptr.get_invert_rtol()

    @property
    def invert_max_iters(self) -> int:
        """Hard iteration cap (termination safeguard) for the inversion."""
        return self._vinet_ptr.get_invert_max_iters()


cdef class InterpolatedEOS(MaterialEOSBase):
    """density(radius) lookup table (e.g. PREM-style profiles)."""

    def __cinit__(self, *args, **kwargs):
        self._interp_ptr = NULL

    def __init__(
            self,
            radius,
            density,
            shear_modulus=None,
            bulk_modulus=None,
            shear_viscosity=None,
            bulk_viscosity=None,
            double thermal_expansion=0.0,
            reference_temperature=None,
            **material):
        cdef c_MaterialEOSConfig config
        cy_fill_material_config(config, material)
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        config.radius      = <vector[double]>radius
        config.density = <vector[double]>density
        if config.radius.size() != config.density.size():
            raise ValueError("radius and density must have the same length.")
        if shear_modulus is not None:
            config.shear_modulus = <vector[double]>shear_modulus
            if config.shear_modulus.size() != config.radius.size():
                raise ValueError("shear_modulus must match radius in length.")
        if bulk_modulus is not None:
            config.bulk_modulus = <vector[double]>bulk_modulus
            if config.bulk_modulus.size() != config.radius.size():
                raise ValueError("bulk_modulus must match radius in length.")
        if shear_viscosity is not None:
            config.shear_viscosity = <vector[double]>shear_viscosity
            if config.shear_viscosity.size() != config.radius.size():
                raise ValueError("shear_viscosity must match radius in length.")
        if bulk_viscosity is not None:
            config.bulk_viscosity = <vector[double]>bulk_viscosity
            if config.bulk_viscosity.size() != config.radius.size():
                raise ValueError("bulk_viscosity must match radius in length.")
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.Interpolated, config)
        self._interp_ptr = <c_InterpolatedEOS*>ptr.get()
        self._eos_ptr    = move(ptr)
        self._ptr        = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._interp_ptr = NULL

    @property
    def num_points(self) -> int:
        """Number of (radius, density) table points."""
        return self._interp_ptr.get_num_points()


# =====================================================================================================================
# Factory
# =====================================================================================================================
# Every config key some material EOS model reads; make_material_eos rejects anything else.
MATERIAL_EOS_CONFIG_KEYS = frozenset({
    "reference_density_kg_m3", "reference_bulk_modulus_pa", "bulk_modulus_derivative", "invert_rtol",
    "invert_max_iters", "radius_m", "density_kg_m3", "shear_modulus_pa", "bulk_modulus_pa",
    "shear_viscosity_pas", "bulk_viscosity_pas", "thermal_expansion_1_k", "reference_temperature_k",
    "shear_modulus_static_pa", "bulk_modulus_static_pa", "shear_viscosity_static_pas", "bulk_viscosity_static_pas",
    "shear_modulus_pressure_derivative", "shear_modulus_temperature_derivative_pa_k",
    "shear_modulus_reference_temperature_k", "thermal_conductivity_w_mk", "heat_capacity_j_kgk",
    # Nested model tables ({"model": ..., ...}), built with make_viscosity / make_partial_melt and attached.
    "shear_viscosity", "bulk_viscosity", "partial_melt"})


def _same_model(str table_name, str model_name) -> bool:
    """Whether two model names, aliases included, resolve to one model; ValueError for a name not in the family."""
    return c_material_eos_model_from_name(table_name.lower().encode("utf-8")) == c_material_eos_model_from_name(model_name.lower().encode("utf-8"))


def make_material_eos(str model_name, dict config=None) -> MaterialEOSBase:
    """Build a material EOS model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        ``"constant"``, ``"bm"``/``"birch_murnaghan"``, ``"vinet"``, or ``"interpolate"`` (case-insensitive).
    config : dict, optional
        ``reference_density_kg_m3``, ``reference_bulk_modulus_pa``, ``bulk_modulus_derivative``, ``invert_rtol``,
        ``invert_max_iters`` (analytic models); ``radius_m`` and ``density_kg_m3`` sequences plus the optional
        ``shear_modulus_pa``, ``bulk_modulus_pa``, ``shear_viscosity_pas``, ``bulk_viscosity_pas`` tables
        (interpolated model); ``thermal_expansion_1_k`` and ``reference_temperature_k`` (every model). The
        material keys, also for every model: ``shear_modulus_static_pa``, ``bulk_modulus_static_pa``,
        ``shear_viscosity_static_pas``, ``bulk_viscosity_static_pas``, the shear law's
        ``shear_modulus_pressure_derivative``, ``shear_modulus_temperature_derivative_pa_k`` and
        ``shear_modulus_reference_temperature_k``, the thermal constants ``thermal_conductivity_w_mk`` and
        ``heat_capacity_j_kgk``, and the nested model tables ``shear_viscosity``,
        ``bulk_viscosity`` and ``partial_melt`` (each a dict with a ``model`` key and that model's own keys).

    Returns
    -------
    MaterialEOSBase
        The concrete model subclass.

    Raises
    ------
    ValueError
        Unknown model name, or a ``config`` key that no material EOS model reads.
    """
    if config is None:
        # No config at all: the defaults of the world-attached path ([layers.default] or [tides] of config_x).
        config = factory_defaults("material", MATERIAL_EOS_CONFIG_KEYS, model_name, _same_model)
    check_config_keys(config, MATERIAL_EOS_CONFIG_KEYS, "material EOS")
    if config is None:
        config = {}
    # Only the supplied keys override the C++ struct defaults.
    cdef c_MaterialEOSConfig cfg
    if "reference_density_kg_m3" in config:
        cfg.reference_density = config["reference_density_kg_m3"]
    if "reference_bulk_modulus_pa" in config:
        cfg.reference_bulk_modulus = config["reference_bulk_modulus_pa"]
    if "bulk_modulus_derivative" in config:
        cfg.bulk_modulus_derivative = config["bulk_modulus_derivative"]
    if "invert_rtol" in config:
        cfg.invert_rtol = config["invert_rtol"]
    if "invert_max_iters" in config:
        cfg.invert_max_iters = config["invert_max_iters"]
    if "thermal_expansion_1_k" in config:
        cfg.thermal_expansion = config["thermal_expansion_1_k"]
    if "reference_temperature_k" in config:
        cfg.reference_temperature = config["reference_temperature_k"]
    if "radius_m" in config:
        cfg.radius = <vector[double]>config["radius_m"]
    if "density_kg_m3" in config:
        cfg.density = <vector[double]>config["density_kg_m3"]
    if "shear_modulus_pa" in config:
        cfg.shear_modulus = <vector[double]>config["shear_modulus_pa"]
    if "bulk_modulus_pa" in config:
        cfg.bulk_modulus = <vector[double]>config["bulk_modulus_pa"]
    if "shear_viscosity_pas" in config:
        cfg.shear_viscosity = <vector[double]>config["shear_viscosity_pas"]
    if "bulk_viscosity_pas" in config:
        cfg.bulk_viscosity = <vector[double]>config["bulk_viscosity_pas"]
    if "shear_modulus_static_pa" in config:
        cfg.shear_modulus_static = config["shear_modulus_static_pa"]
    if "bulk_modulus_static_pa" in config:
        cfg.bulk_modulus_static = config["bulk_modulus_static_pa"]
    if "shear_viscosity_static_pas" in config:
        cfg.shear_viscosity_static = config["shear_viscosity_static_pas"]
    if "bulk_viscosity_static_pas" in config:
        cfg.bulk_viscosity_static = config["bulk_viscosity_static_pas"]
    if "shear_modulus_pressure_derivative" in config:
        cfg.shear_modulus_pressure_derivative = config["shear_modulus_pressure_derivative"]
    if "shear_modulus_temperature_derivative_pa_k" in config:
        cfg.shear_modulus_temperature_derivative = config["shear_modulus_temperature_derivative_pa_k"]
    if "shear_modulus_reference_temperature_k" in config:
        cfg.shear_modulus_reference_temperature = config["shear_modulus_reference_temperature_k"]
    if "thermal_conductivity_w_mk" in config:
        cfg.thermal_conductivity = config["thermal_conductivity_w_mk"]
    if "heat_capacity_j_kgk" in config:
        cfg.heat_capacity = config["heat_capacity_j_kgk"]

    cdef c_MaterialEOSModel model = c_material_eos_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(model, cfg)

    cdef ConstantDensityEOS constant_eos
    cdef BirchMurnaghanEOS  bm_eos
    cdef VinetEOS           vinet_eos
    cdef InterpolatedEOS    interp_eos
    if model == c_MaterialEOSModel.Constant:
        constant_eos = ConstantDensityEOS.__new__(ConstantDensityEOS)
        constant_eos._constant_ptr = <c_ConstantDensityEOS*>ptr.get()
        constant_eos._eos_ptr = move(ptr)
        constant_eos._ptr = <c_TidalPyBaseClass*>constant_eos._eos_ptr.get()
        return _attach_material_models(constant_eos, config)
    elif model == c_MaterialEOSModel.BirchMurnaghan:
        bm_eos = BirchMurnaghanEOS.__new__(BirchMurnaghanEOS)
        bm_eos._bm_ptr = <c_BirchMurnaghanEOS*>ptr.get()
        bm_eos._eos_ptr = move(ptr)
        bm_eos._ptr = <c_TidalPyBaseClass*>bm_eos._eos_ptr.get()
        return _attach_material_models(bm_eos, config)
    elif model == c_MaterialEOSModel.Vinet:
        vinet_eos = VinetEOS.__new__(VinetEOS)
        vinet_eos._vinet_ptr = <c_VinetEOS*>ptr.get()
        vinet_eos._eos_ptr = move(ptr)
        vinet_eos._ptr = <c_TidalPyBaseClass*>vinet_eos._eos_ptr.get()
        return _attach_material_models(vinet_eos, config)
    else:
        interp_eos = InterpolatedEOS.__new__(InterpolatedEOS)
        interp_eos._interp_ptr = <c_InterpolatedEOS*>ptr.get()
        interp_eos._eos_ptr = move(ptr)
        interp_eos._ptr = <c_TidalPyBaseClass*>interp_eos._eos_ptr.get()
        return _attach_material_models(interp_eos, config)


def _attach_material_models(MaterialEOSBase eos, dict config) -> MaterialEOSBase:
    """Build and attach the nested viscosity and partial-melt tables of a material config."""
    # Imported here: these modules sit beside this one and nothing at import time needs them.
    from TidalPy.viscosity_x.viscosity import make_viscosity
    from TidalPy.partial_melt_x.partial_melt import make_partial_melt

    cdef str key
    cdef object setter
    cdef object maker
    cdef dict section
    cdef str model_name
    # `model` is whichever model class the maker returns (a viscosity or a partial-melt model), so it stays
    # `object`; the two makers do not share a base class.
    cdef object model
    for key, setter, maker in (
            ("shear_viscosity", eos.set_shear_viscosity, make_viscosity),
            ("bulk_viscosity", eos.set_bulk_viscosity, make_viscosity),
            ("partial_melt", eos.set_partial_melt, make_partial_melt)):
        section = config.get(key)
        if section is None:
            continue
        section = dict(section)
        model_name = section.pop("model", None)
        if model_name is None:
            raise ValueError(f"The material's {key!r} table needs a 'model' key.")
        try:
            model = maker(model_name, section)
        except ValueError as error:
            # Name the nested table so a rejected key can be found in the source file.
            raise ValueError(f"[{key}] {error}") from error
        setter(model)
    return eos
