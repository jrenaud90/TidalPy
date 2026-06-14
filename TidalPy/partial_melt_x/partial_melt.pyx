# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
partial_melt.pyx
Cython/Python wrapper for TidalPy's partial-melt model hierarchy.

A partial-melt model maps a material's pre-melt (solid) viscosity and shear
modulus, plus its temperature, to the post-melt viscosity and shear modulus (melt
weakening), and reports the volumetric melt fraction. These quantities are
frequency-independent and feed the downstream rheology (complex modulus) step.

Models:
- OffPartialMelt     (alias "none")    — no weakening (returns pre-melt).
- SpohnPartialMelt   (alias "fischer") — Fischer & Spohn (1990) temperature law.
- HenningPartialMelt                   — Henning (2009/2010) three-regime law.
"""

from libcpp.string cimport string
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
# PartialMeltBase
# =====================================================================================================================

cdef class PartialMeltBase(PhysicsBase):
    """Abstract base for partial-melt models. Instantiate a concrete subclass."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_PartialMeltBase> auto-inits to nullptr

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "PartialMeltBase is abstract; instantiate a concrete model "
            "(OffPartialMelt, SpohnPartialMelt, HenningPartialMelt).")

    def __dealloc__(self):
        self._melt_ptr.reset()
        self._ptr = NULL

    @property
    def solidus(self) -> float:
        """Solidus temperature [K]."""
        return self._melt_ptr.get().get_solidus()

    @property
    def liquidus(self) -> float:
        """Liquidus temperature [K]."""
        return self._melt_ptr.get().get_liquidus()

    @property
    def liquid_shear(self) -> float:
        """Shear modulus of the fully-molten material [Pa]."""
        return self._melt_ptr.get().get_liquid_shear()

    def calc_melt_fraction(self, double temperature_k) -> float:
        """Volumetric melt fraction phi in [0, 1] from temperature [K]."""
        return self._melt_ptr.get().calc_melt_fraction(temperature_k)

    def calc_partial_melt(
            self,
            double temperature_k,
            double premelt_viscosity,
            double premelt_shear,
            double liquid_viscosity) -> tuple:
        """Post-melt strength from the pre-melt state (all MKS).

        Returns
        -------
        (melt_fraction, postmelt_viscosity, postmelt_shear_modulus) : tuple of float
        """
        cdef c_PartialMeltInputs inputs
        inputs.temperature_k     = temperature_k
        inputs.premelt_viscosity = premelt_viscosity
        inputs.premelt_shear     = premelt_shear
        inputs.liquid_viscosity  = liquid_viscosity
        cdef c_PartialMeltResult result = self._melt_ptr.get().calc_partial_melt(inputs)
        return (result.melt_fraction, result.postmelt_viscosity, result.postmelt_shear_modulus)

    cpdef dict get_config_dict(self):
        """Return the model name as a config dict (subclasses add parameters)."""
        cdef c_PartialMeltBase* p = self._melt_ptr.get()
        return {
            "model":           self.model_name,
            "solidus_k":       p.get_solidus(),
            "liquidus_k":      p.get_liquidus(),
            "liquid_shear_pa": p.get_liquid_shear(),
        }


# =====================================================================================================================
# Partial-melt models
# =====================================================================================================================
cdef class OffPartialMelt(PartialMeltBase):
    """No melt weakening; post-melt strength equals pre-melt."""

    def __cinit__(self, *args, **kwargs):
        self._off_ptr = NULL

    def __init__(self, double solidus_k=1600.0, double liquidus_k=2000.0,
                 double liquid_shear_pa=1.0e-5):
        cdef c_PartialMeltConfig config
        config.solidus_k       = solidus_k
        config.liquidus_k      = liquidus_k
        config.liquid_shear_pa = liquid_shear_pa
        cdef unique_ptr[c_PartialMeltBase] ptr = c_find_partial_melt(
            c_PartialMeltModel.Off, config)
        self._off_ptr  = <c_OffPartialMelt*>ptr.get()
        self._melt_ptr = move(ptr)
        self._ptr      = <c_TidalPyBaseClass*>self._melt_ptr.get()

    def __dealloc__(self):
        self._off_ptr = NULL


cdef class SpohnPartialMelt(PartialMeltBase):
    """Fischer & Spohn (1990) temperature-based melt law."""

    def __cinit__(self, *args, **kwargs):
        self._spohn_ptr = NULL

    def __init__(
            self,
            double solidus_k=1600.0,
            double liquidus_k=2000.0,
            double liquid_shear_pa=1.0e-5,
            double fs_visc_power_slope=27000.0,
            double fs_visc_power_phase=1.0,
            double fs_shear_power_slope=82000.0,
            double fs_shear_power_phase=40.6):
        cdef c_PartialMeltConfig config
        config.solidus_k            = solidus_k
        config.liquidus_k           = liquidus_k
        config.liquid_shear_pa      = liquid_shear_pa
        config.fs_visc_power_slope  = fs_visc_power_slope
        config.fs_visc_power_phase  = fs_visc_power_phase
        config.fs_shear_power_slope = fs_shear_power_slope
        config.fs_shear_power_phase = fs_shear_power_phase
        cdef unique_ptr[c_PartialMeltBase] ptr = c_find_partial_melt(c_PartialMeltModel.Spohn, config)
        self._spohn_ptr = <c_SpohnPartialMelt*>ptr.get()
        self._melt_ptr  = move(ptr)
        self._ptr       = <c_TidalPyBaseClass*>self._melt_ptr.get()

    def __dealloc__(self):
        self._spohn_ptr = NULL

    cpdef dict get_config_dict(self):
        d = PartialMeltBase.get_config_dict(self)
        d["fs_visc_power_slope"]  = self._spohn_ptr.get_visc_power_slope()
        d["fs_visc_power_phase"]  = self._spohn_ptr.get_visc_power_phase()
        d["fs_shear_power_slope"] = self._spohn_ptr.get_shear_power_slope()
        d["fs_shear_power_phase"] = self._spohn_ptr.get_shear_power_phase()
        return d


cdef class HenningPartialMelt(PartialMeltBase):
    """Henning (2009/2010) three-regime melt weakening."""

    def __cinit__(self, *args, **kwargs):
        self._henning_ptr = NULL

    def __init__(
            self,
            double solidus_k=1600.0,
            double liquidus_k=2000.0,
            double liquid_shear_pa=1.0e-5,
            double crit_melt_frac=0.5,
            double crit_melt_frac_width=0.05,
            double hn_visc_slope_1=13.5,
            double hn_visc_falloff_slope=370.0,
            double hn_shear_param_1=40000.0,
            double hn_shear_param_2=25.0,
            double hn_shear_falloff_slope=700.0):
        cdef c_PartialMeltConfig config
        config.solidus_k              = solidus_k
        config.liquidus_k             = liquidus_k
        config.liquid_shear_pa        = liquid_shear_pa
        config.crit_melt_frac         = crit_melt_frac
        config.crit_melt_frac_width   = crit_melt_frac_width
        config.hn_visc_slope_1        = hn_visc_slope_1
        config.hn_visc_falloff_slope  = hn_visc_falloff_slope
        config.hn_shear_param_1       = hn_shear_param_1
        config.hn_shear_param_2       = hn_shear_param_2
        config.hn_shear_falloff_slope = hn_shear_falloff_slope
        cdef unique_ptr[c_PartialMeltBase] ptr = c_find_partial_melt(c_PartialMeltModel.Henning, config)
        self._henning_ptr = <c_HenningPartialMelt*>ptr.get()
        self._melt_ptr    = move(ptr)
        self._ptr         = <c_TidalPyBaseClass*>self._melt_ptr.get()

    def __dealloc__(self):
        self._henning_ptr = NULL

    cpdef dict get_config_dict(self):
        d = PartialMeltBase.get_config_dict(self)
        d["crit_melt_frac"]         = self._henning_ptr.get_crit_melt_frac()
        d["crit_melt_frac_width"]   = self._henning_ptr.get_crit_melt_frac_width()
        d["hn_visc_slope_1"]        = self._henning_ptr.get_visc_slope_1()
        d["hn_visc_falloff_slope"]  = self._henning_ptr.get_visc_falloff_slope()
        d["hn_shear_param_1"]       = self._henning_ptr.get_shear_param_1()
        d["hn_shear_param_2"]       = self._henning_ptr.get_shear_param_2()
        d["hn_shear_falloff_slope"] = self._henning_ptr.get_shear_falloff_slope()
        return d


# =====================================================================================================================
# Factory
# =====================================================================================================================
def make_partial_melt(str model_name, dict config=None) -> PartialMeltBase:
    """Build a partial-melt model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        One of ``"off"``/``"none"``, ``"spohn"``/``"fischer"``, ``"henning"``
        (case-insensitive; aliases accepted).
    config : dict, optional
        Model parameters (solidus_k, liquidus_k, liquid_shear_pa, plus the
        Spohn/Henning scalars). Absent keys fall back to the C++ defaults.

    Returns
    -------
    PartialMeltBase
        The concrete model subclass.

    Raises
    ------
    ValueError
        If the model name is unknown.
    """
    if config is None:
        config = {}
    # A default-constructed config carries the C++ defaults; only override the
    # fields the caller actually supplies (single source of truth: the C++ struct).
    cdef c_PartialMeltConfig cfg
    if "solidus_k" in config:              cfg.solidus_k              = config["solidus_k"]
    if "liquidus_k" in config:             cfg.liquidus_k             = config["liquidus_k"]
    if "liquid_shear_pa" in config:        cfg.liquid_shear_pa        = config["liquid_shear_pa"]
    if "fs_visc_power_slope" in config:    cfg.fs_visc_power_slope    = config["fs_visc_power_slope"]
    if "fs_visc_power_phase" in config:    cfg.fs_visc_power_phase    = config["fs_visc_power_phase"]
    if "fs_shear_power_slope" in config:   cfg.fs_shear_power_slope   = config["fs_shear_power_slope"]
    if "fs_shear_power_phase" in config:   cfg.fs_shear_power_phase   = config["fs_shear_power_phase"]
    if "crit_melt_frac" in config:         cfg.crit_melt_frac         = config["crit_melt_frac"]
    if "crit_melt_frac_width" in config:   cfg.crit_melt_frac_width   = config["crit_melt_frac_width"]
    if "hn_visc_slope_1" in config:        cfg.hn_visc_slope_1        = config["hn_visc_slope_1"]
    if "hn_visc_falloff_slope" in config:  cfg.hn_visc_falloff_slope  = config["hn_visc_falloff_slope"]
    if "hn_shear_param_1" in config:       cfg.hn_shear_param_1       = config["hn_shear_param_1"]
    if "hn_shear_param_2" in config:       cfg.hn_shear_param_2       = config["hn_shear_param_2"]
    if "hn_shear_falloff_slope" in config: cfg.hn_shear_falloff_slope = config["hn_shear_falloff_slope"]

    cdef c_PartialMeltModel model = c_partial_melt_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_PartialMeltBase] ptr = c_find_partial_melt(model, cfg)

    cdef OffPartialMelt     off_melt
    cdef SpohnPartialMelt   spohn_melt
    cdef HenningPartialMelt henning_melt
    if model == c_PartialMeltModel.Off:
        off_melt = OffPartialMelt.__new__(OffPartialMelt)
        off_melt._off_ptr  = <c_OffPartialMelt*>ptr.get()
        off_melt._melt_ptr = move(ptr)
        off_melt._ptr      = <c_TidalPyBaseClass*>off_melt._melt_ptr.get()
        return off_melt
    elif model == c_PartialMeltModel.Spohn:
        spohn_melt = SpohnPartialMelt.__new__(SpohnPartialMelt)
        spohn_melt._spohn_ptr = <c_SpohnPartialMelt*>ptr.get()
        spohn_melt._melt_ptr  = move(ptr)
        spohn_melt._ptr       = <c_TidalPyBaseClass*>spohn_melt._melt_ptr.get()
        return spohn_melt
    else:
        henning_melt = HenningPartialMelt.__new__(HenningPartialMelt)
        henning_melt._henning_ptr = <c_HenningPartialMelt*>ptr.get()
        henning_melt._melt_ptr    = move(ptr)
        henning_melt._ptr         = <c_TidalPyBaseClass*>henning_melt._melt_ptr.get()
        return henning_melt
