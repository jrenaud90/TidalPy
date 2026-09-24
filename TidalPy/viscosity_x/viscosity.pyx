# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's viscosity models."""

from libcpp.string cimport string
from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass, cy_resolve_factory_config

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class ViscosityBase(PhysicsBase):
    """Abstract base for viscosity models. Instantiate a concrete subclass."""

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "ViscosityBase is abstract; instantiate a concrete model "
            "(ArrheniusViscosity, ReferenceViscosity, ConstantViscosity).")

    def __dealloc__(self):
        self._visc_ptr.reset()
        self._ptr = NULL

    cdef void _adopt(self, unique_ptr[c_ViscosityBase]& model) noexcept:
        """Take ownership of ``model``; the inherited ``_ptr`` observes it."""
        self._visc_ptr = move(model)
        self._ptr = <c_TidalPyBaseClass*>self._visc_ptr.get()

    def calc_viscosity(self, double temperature, double pressure=0.0) -> float:
        """Dynamic viscosity [Pa s] at the given temperature [K] and pressure [Pa]."""
        self._check_ptr()
        return self._visc_ptr.get().calc_viscosity(temperature, pressure)


cdef class ConstantViscosity(ViscosityBase):
    """Viscosity independent of temperature and pressure."""

    def __init__(self, double reference_viscosity=1.0e22):
        cdef c_ViscosityConfig config
        config.reference_viscosity = reference_viscosity
        cdef unique_ptr[c_ViscosityBase] model = c_find_viscosity(c_ViscosityModel.Constant, config)
        self._adopt(model)

    @property
    def reference_viscosity(self) -> float:
        """Reference (constant) viscosity [Pa s]."""
        self._check_ptr()
        return (<c_ConstantViscosity*>self._visc_ptr.get()).get_reference_viscosity()


cdef class ReferenceViscosity(ViscosityBase):
    """Relative-activation law: eta = eta_ref * exp(((E_a + P*V_a)/R)*(1/T - 1/T_ref))."""

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
        cdef unique_ptr[c_ViscosityBase] model = c_find_viscosity(c_ViscosityModel.Reference, config)
        self._adopt(model)

    @property
    def reference_viscosity(self) -> float:
        """Reference viscosity [Pa s]."""
        self._check_ptr()
        return (<c_ReferenceViscosity*>self._visc_ptr.get()).get_reference_viscosity()

    @property
    def reference_temperature(self) -> float:
        """Reference temperature [K]."""
        self._check_ptr()
        return (<c_ReferenceViscosity*>self._visc_ptr.get()).get_reference_temperature()

    @property
    def molar_activation_energy(self) -> float:
        """Molar activation energy E_a [J/mol]."""
        self._check_ptr()
        return (<c_ReferenceViscosity*>self._visc_ptr.get()).get_molar_activation_energy()

    @property
    def molar_activation_volume(self) -> float:
        """Molar activation volume V_a [m^3/mol]."""
        self._check_ptr()
        return (<c_ReferenceViscosity*>self._visc_ptr.get()).get_molar_activation_volume()


cdef class ArrheniusViscosity(ViscosityBase):
    """Arrhenius flow law: eta = A * sigma^(1-n) * d^m * exp((E_a + P*V_a)/(R*T))."""

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
        cdef unique_ptr[c_ViscosityBase] model = c_find_viscosity(c_ViscosityModel.Arrhenius, config)
        self._adopt(model)

    @property
    def arrhenius_coeff(self) -> float:
        """Pre-exponential coefficient A."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_arrhenius_coeff()

    @property
    def stress(self) -> float:
        """Applied shear stress sigma [Pa]; the stress term drops out when ``stress_expo`` is 1."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_stress()

    @property
    def stress_expo(self) -> float:
        """Stress exponent n: 1 for diffusion creep, above 1 for dislocation creep."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_stress_expo()

    @property
    def grain_size(self) -> float:
        """Grain size d [m]."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_grain_size()

    @property
    def grain_size_expo(self) -> float:
        """Grain-size exponent m; 0 removes the grain-size dependence."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_grain_size_expo()

    @property
    def molar_activation_energy(self) -> float:
        """Molar activation energy E_a [J/mol]."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_molar_activation_energy()

    @property
    def molar_activation_volume(self) -> float:
        """Molar activation volume V_a [m^3/mol]."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_molar_activation_volume()

    @property
    def additional_temp_dependence(self) -> bool:
        """Whether the law is multiplied by an additional factor of T."""
        self._check_ptr()
        return (<c_ArrheniusViscosity*>self._visc_ptr.get()).get_additional_temp_dependence()


# Every config key any viscosity model reads; make_viscosity rejects anything else.
VISCOSITY_CONFIG_KEYS = frozenset({
    "reference_viscosity_pas", "reference_temperature_k", "molar_activation_energy_j_mol",
    "molar_activation_volume_m3_mol", "arrhenius_coeff", "stress_pa", "stress_expo", "grain_size_m",
    "grain_size_expo", "additional_temp_dependence"})


# The wrapper class of each c_ViscosityModel, in enum order.
_VISCOSITY_CLASSES = (ArrheniusViscosity, ReferenceViscosity, ConstantViscosity)


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return (
        c_viscosity_model_from_name(table_name.encode("utf-8"))
        == c_viscosity_model_from_name(model_name.encode("utf-8")))


def make_viscosity(str model_name, dict config=None) -> ViscosityBase:
    """Build a viscosity model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        ``"arrhenius"``/``"arr"``, ``"reference"``/``"ref"``, or ``"constant"``/``"const"``.
    config : dict, optional
        Model parameters, keyed with their units (see ``VISCOSITY_CONFIG_KEYS``). Absent keys fall back to
        the C++ defaults.

    Returns
    -------
    ViscosityBase

    Raises
    ------
    ValueError
        Unknown model name, or a config key that no viscosity model reads.
    """
    # None falls back to the same defaults the world-attached path uses.
    config = cy_resolve_factory_config(
        config, "material.shear_viscosity", VISCOSITY_CONFIG_KEYS, model_name, _same_model, "viscosity")
    # The default-constructed config carries the C++ defaults, so only override what the caller gave.
    cdef c_ViscosityConfig cfg
    cfg.reference_viscosity        = config.get("reference_viscosity_pas", cfg.reference_viscosity)
    cfg.reference_temperature      = config.get("reference_temperature_k", cfg.reference_temperature)
    cfg.molar_activation_energy    = config.get("molar_activation_energy_j_mol", cfg.molar_activation_energy)
    cfg.molar_activation_volume    = config.get("molar_activation_volume_m3_mol", cfg.molar_activation_volume)
    cfg.arrhenius_coeff            = config.get("arrhenius_coeff", cfg.arrhenius_coeff)
    cfg.stress                     = config.get("stress_pa", cfg.stress)
    cfg.stress_expo                = config.get("stress_expo", cfg.stress_expo)
    cfg.grain_size                 = config.get("grain_size_m", cfg.grain_size)
    cfg.grain_size_expo            = config.get("grain_size_expo", cfg.grain_size_expo)
    cfg.additional_temp_dependence = bool(config.get("additional_temp_dependence", cfg.additional_temp_dependence))

    cdef c_ViscosityModel model = c_viscosity_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_ViscosityBase] ptr = c_find_viscosity(model, cfg)
    wrapper_class = _VISCOSITY_CLASSES[<int>model]
    cdef ViscosityBase wrapper = wrapper_class.__new__(wrapper_class)
    wrapper._adopt(ptr)
    return wrapper
