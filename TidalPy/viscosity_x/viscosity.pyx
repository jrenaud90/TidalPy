# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
viscosity.pyx
Cython/Python wrapper for TidalPy's viscosity model hierarchy.

A viscosity model returns a material's dynamic viscosity [Pa·s] as a function of
temperature [K] and pressure [Pa] — the pre-melt ("solid") viscosity that the
partial-melt step subsequently weakens.

Models:
- ArrheniusViscosity  (alias "arr")   — Arrhenius flow law.
- ReferenceViscosity  (alias "ref")   — relative-activation law.
- ConstantViscosity   (alias "const") — temperature/pressure independent.
"""

from libcpp.string cimport string
from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# ViscosityBase
# =====================================================================================================================

cdef class ViscosityBase(PhysicsBase):
    """Abstract base for viscosity models. Instantiate a concrete subclass."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_ViscosityBase> auto-inits to nullptr

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "ViscosityBase is abstract; instantiate a concrete model "
            "(ArrheniusViscosity, ReferenceViscosity, ConstantViscosity).")

    def __dealloc__(self):
        self._visc_ptr.reset()
        self._ptr = NULL

    def calc_viscosity(self, double temperature_k, double pressure_pa=0.0) -> float:
        """Dynamic viscosity [Pa·s] at the given temperature [K] and pressure [Pa]."""
        return self._visc_ptr.get().calc_viscosity(temperature_k, pressure_pa)

    cpdef dict get_config_dict(self):
        """Return the model name as a config dict (subclasses add parameters)."""
        return {"model": self.model_name}


# =====================================================================================================================
# Viscosity models
# =====================================================================================================================
cdef class ConstantViscosity(ViscosityBase):
    """Viscosity independent of temperature and pressure."""

    def __cinit__(self, *args, **kwargs):
        self._constant_ptr = NULL

    def __init__(self, double reference_viscosity=1.0e22):
        cdef c_ViscosityConfig config
        config.reference_viscosity = reference_viscosity
        cdef unique_ptr[c_ViscosityBase] ptr = c_find_viscosity(c_ViscosityModel.Constant, config)
        self._constant_ptr = <c_ConstantViscosity*>ptr.get()
        self._visc_ptr     = move(ptr)
        self._ptr          = <c_TidalPyBaseClass*>self._visc_ptr.get()

    def __dealloc__(self):
        self._constant_ptr = NULL

    @property
    def reference_viscosity(self) -> float:
        """Reference (constant) viscosity [Pa·s]."""
        return self._constant_ptr.get_reference_viscosity()

    cpdef dict get_config_dict(self):
        d = ViscosityBase.get_config_dict(self)
        d["reference_viscosity"] = self._constant_ptr.get_reference_viscosity()
        return d


cdef class ReferenceViscosity(ViscosityBase):
    """Relative-activation law: eta = eta_ref * exp(((E_a + P*V_a)/R)*(1/T - 1/T_ref))."""

    def __cinit__(self, *args, **kwargs):
        self._ref_ptr = NULL

    def __init__(
            self,
            double reference_viscosity=1.0e22,
            double reference_temperature=1000.0,
            double molar_activation_energy=3.0e5,
            double molar_activation_volume=0.0):
        cdef c_ViscosityConfig config
        config.reference_viscosity     = reference_viscosity
        config.reference_temperature   = reference_temperature
        config.molar_activation_energy = molar_activation_energy
        config.molar_activation_volume = molar_activation_volume
        cdef unique_ptr[c_ViscosityBase] ptr = c_find_viscosity(c_ViscosityModel.Reference, config)
        self._ref_ptr  = <c_ReferenceViscosity*>ptr.get()
        self._visc_ptr = move(ptr)
        self._ptr      = <c_TidalPyBaseClass*>self._visc_ptr.get()

    def __dealloc__(self):
        self._ref_ptr = NULL

    @property
    def reference_viscosity(self) -> float:
        """Reference viscosity [Pa·s]."""
        return self._ref_ptr.get_reference_viscosity()

    @property
    def reference_temperature(self) -> float:
        """Reference temperature [K]."""
        return self._ref_ptr.get_reference_temperature()

    @property
    def molar_activation_energy(self) -> float:
        """Molar activation energy E_a [J/mol]."""
        return self._ref_ptr.get_molar_activation_energy()

    @property
    def molar_activation_volume(self) -> float:
        """Molar activation volume V_a [m^3/mol]."""
        return self._ref_ptr.get_molar_activation_volume()

    cpdef dict get_config_dict(self):
        d = ViscosityBase.get_config_dict(self)
        d["reference_viscosity"]     = self._ref_ptr.get_reference_viscosity()
        d["reference_temperature"]   = self._ref_ptr.get_reference_temperature()
        d["molar_activation_energy"] = self._ref_ptr.get_molar_activation_energy()
        d["molar_activation_volume"] = self._ref_ptr.get_molar_activation_volume()
        return d


cdef class ArrheniusViscosity(ViscosityBase):
    """Arrhenius flow law: eta = A * sigma^(1-n) * d^m * exp((E_a + P*V_a)/(R*T))."""

    def __cinit__(self, *args, **kwargs):
        self._arr_ptr = NULL

    def __init__(
            self,
            double arrhenius_coeff=1.0,
            double stress=1.0,
            double stress_expo=1.0,
            double grain_size=1.0e-3,
            double grain_size_expo=0.0,
            double molar_activation_energy=3.0e5,
            double molar_activation_volume=0.0,
            cpp_bool additional_temp_dependence=False):
        cdef c_ViscosityConfig config
        config.arrhenius_coeff            = arrhenius_coeff
        config.stress                     = stress
        config.stress_expo                = stress_expo
        config.grain_size                 = grain_size
        config.grain_size_expo            = grain_size_expo
        config.molar_activation_energy    = molar_activation_energy
        config.molar_activation_volume    = molar_activation_volume
        config.additional_temp_dependence = additional_temp_dependence
        cdef unique_ptr[c_ViscosityBase] ptr = c_find_viscosity(c_ViscosityModel.Arrhenius, config)
        self._arr_ptr  = <c_ArrheniusViscosity*>ptr.get()
        self._visc_ptr = move(ptr)
        self._ptr      = <c_TidalPyBaseClass*>self._visc_ptr.get()

    def __dealloc__(self):
        self._arr_ptr = NULL

    @property
    def arrhenius_coeff(self) -> float:
        """Pre-exponential coefficient A."""
        return self._arr_ptr.get_arrhenius_coeff()

    @property
    def molar_activation_energy(self) -> float:
        """Molar activation energy E_a [J/mol]."""
        return self._arr_ptr.get_molar_activation_energy()

    @property
    def additional_temp_dependence(self) -> bool:
        """Whether the law is multiplied by an additional factor of T."""
        return self._arr_ptr.get_additional_temp_dependence()

    cpdef dict get_config_dict(self):
        d = ViscosityBase.get_config_dict(self)
        d["arrhenius_coeff"]            = self._arr_ptr.get_arrhenius_coeff()
        d["stress"]                     = self._arr_ptr.get_stress()
        d["stress_expo"]                = self._arr_ptr.get_stress_expo()
        d["grain_size"]                 = self._arr_ptr.get_grain_size()
        d["grain_size_expo"]            = self._arr_ptr.get_grain_size_expo()
        d["molar_activation_energy"]    = self._arr_ptr.get_molar_activation_energy()
        d["molar_activation_volume"]    = self._arr_ptr.get_molar_activation_volume()
        d["additional_temp_dependence"] = bool(self._arr_ptr.get_additional_temp_dependence())
        return d


# =====================================================================================================================
# Factory
# =====================================================================================================================
def make_viscosity(str model_name, dict config=None) -> ViscosityBase:
    """Build a viscosity model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        One of ``"arrhenius"``/``"arr"``, ``"reference"``/``"ref"``,
        ``"constant"``/``"const"`` (case-insensitive; aliases accepted).
    config : dict, optional
        Model parameters; absent keys fall back to the C++ defaults.

    Returns
    -------
    ViscosityBase
        The concrete model subclass.

    Raises
    ------
    ValueError
        If the model name is unknown.
    """
    if config is None:
        config = {}
    cdef c_ViscosityConfig cfg
    if "reference_viscosity" in config:     cfg.reference_viscosity     = config["reference_viscosity"]
    if "reference_temperature" in config:   cfg.reference_temperature   = config["reference_temperature"]
    if "molar_activation_energy" in config: cfg.molar_activation_energy = config["molar_activation_energy"]
    if "molar_activation_volume" in config: cfg.molar_activation_volume = config["molar_activation_volume"]
    if "arrhenius_coeff" in config:         cfg.arrhenius_coeff         = config["arrhenius_coeff"]
    if "stress" in config:                  cfg.stress                  = config["stress"]
    if "stress_expo" in config:             cfg.stress_expo             = config["stress_expo"]
    if "grain_size" in config:              cfg.grain_size              = config["grain_size"]
    if "grain_size_expo" in config:         cfg.grain_size_expo         = config["grain_size_expo"]
    if "additional_temp_dependence" in config:
        cfg.additional_temp_dependence = bool(config["additional_temp_dependence"])

    cdef c_ViscosityModel model = c_viscosity_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_ViscosityBase] ptr = c_find_viscosity(model, cfg)

    cdef ArrheniusViscosity arr_visc
    cdef ReferenceViscosity ref_visc
    cdef ConstantViscosity  const_visc
    if model == c_ViscosityModel.Arrhenius:
        arr_visc = ArrheniusViscosity.__new__(ArrheniusViscosity)
        arr_visc._arr_ptr  = <c_ArrheniusViscosity*>ptr.get()
        arr_visc._visc_ptr = move(ptr)
        arr_visc._ptr      = <c_TidalPyBaseClass*>arr_visc._visc_ptr.get()
        return arr_visc
    elif model == c_ViscosityModel.Reference:
        ref_visc = ReferenceViscosity.__new__(ReferenceViscosity)
        ref_visc._ref_ptr  = <c_ReferenceViscosity*>ptr.get()
        ref_visc._visc_ptr = move(ptr)
        ref_visc._ptr      = <c_TidalPyBaseClass*>ref_visc._visc_ptr.get()
        return ref_visc
    else:
        const_visc = ConstantViscosity.__new__(ConstantViscosity)
        const_visc._constant_ptr = <c_ConstantViscosity*>ptr.get()
        const_visc._visc_ptr     = move(ptr)
        const_visc._ptr          = <c_TidalPyBaseClass*>const_visc._visc_ptr.get()
        return const_visc
