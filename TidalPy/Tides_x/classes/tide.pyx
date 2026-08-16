# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
tide.pyx
Cython/Python wrapper for TidalPy's global (1D) tidal dissipation model hierarchy.

A tide model converts a per-mode Love number into the dissipation multiplier -Im[k_l] 
used by the global mode collapse. The four models:

- RheologyTide  (alias "rheology")              — k_l from the radial solver.
- FixedQTide    (alias "cpl"/"fixed_q")         — constant phase lag, k_l*(1 - i/Q_l).
- FixedLagTide  (alias "ctl"/"fixed_dt")        — constant time lag,  k_l*(1 - i*w*dt_l).
- CTLQTide      (alias "ctl_q"/"fixed_dt_q")    — k_l*(1 - i*w*dt_l/Q_l).

Per-degree fixed parameters (fixed_k, fixed_q, fixed_dt) are supplied as lists indexed
from degree l = 2 (index 0 -> l=2). Supported degrees are l = 2..10.
"""

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Supported degree range (matches the C++ tide_.hpp slots; l = 2..10).
cdef int C_TIDE_MIN_DEGREE = 2
cdef int C_TIDE_MAX_DEGREE = 10


# =====================================================================================================================
# Config helpers
# =====================================================================================================================
cdef vector[double] _to_double_vector(object values) except *:
    """Convert a Python iterable of floats (indexed from l=2) into a std::vector[double]."""
    cdef vector[double] out
    if values is None:
        return out
    for value in values:
        out.push_back(<double>value)
    return out


cdef c_TideModelConfig _build_tide_config(dict config) except *:
    """Build a c_TideModelConfig from a config dict with optional per-degree list keys."""
    cdef c_TideModelConfig cfg
    if config is None:
        return cfg
    if "fixed_k" in config:
        cfg.fixed_k = _to_double_vector(config["fixed_k"])
    if "fixed_q" in config:
        cfg.fixed_q = _to_double_vector(config["fixed_q"])
    if "fixed_dt" in config:
        cfg.fixed_dt = _to_double_vector(config["fixed_dt"])
    return cfg


cdef list _per_degree_list(object getter):
    """Build a per-degree list [l=2 .. l=10] from a C++ get_fixed_*(degree_l) callable."""
    cdef list out = []
    cdef int degree_l
    for degree_l in range(C_TIDE_MIN_DEGREE, C_TIDE_MAX_DEGREE + 1):
        out.append(getter(degree_l))
    return out


# =====================================================================================================================
# TideBase
# =====================================================================================================================
cdef class TideBase(PhysicsBase):
    """Abstract base for global tide models. Instantiate a concrete subclass."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_TideBase> auto-inits to nullptr

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "TideBase is abstract; instantiate a concrete model "
            "(RheologyTide, FixedQTide, FixedLagTide, CTLQTide).")

    def __dealloc__(self):
        self._tide_ptr.reset()
        self._ptr = NULL

    def calc_love_numbers(self, int degree_l, double frequency, LoveNumbers solver_love=None) -> LoveNumbers:
        """Full complex Love-number suite (k, h, l) at the tidal frequency [rad s-1].

        For analytic models h and l are NaN (no radial solution); the rheology model
        returns the supplied ``solver_love`` (a ``LoveNumbers``) unchanged. ``solver_love``
        is ignored by the analytic models and defaults to zeros when omitted.
        """
        cdef c_LoveNumbers solver_c
        if solver_love is not None:
            solver_c = solver_love._love
        cdef c_LoveNumbers result = self._tide_ptr.get().calc_love_numbers(degree_l, frequency, solver_c)
        cdef LoveNumbers out = LoveNumbers.__new__(LoveNumbers)
        out._love = result
        return out

    def calc_neg_imk(self, int degree_l, double frequency, LoveNumbers solver_love=None) -> float:
        """-Im[k_l] at the tidal frequency [rad s-1] (the mode-collapse dissipation multiplier)."""
        cdef c_LoveNumbers solver_c
        if solver_love is not None:
            solver_c = solver_love._love
        return self._tide_ptr.get().calc_neg_imk(degree_l, frequency, solver_c)

    @property
    def needs_radial_solve(self) -> bool:
        """Whether this model requires the radial solver to supply k_l (rheology)."""
        return bool(self._tide_ptr.get().needs_radial_solve())

    cpdef dict get_config_dict(self):
        """Return the model name as a config dict (subclasses add parameters)."""
        return {"model": self.model_name}


# =====================================================================================================================
# Tide models
# =====================================================================================================================
cdef class RheologyTide(TideBase):
    """Rheology-based tide: k_l comes from the radial solver (frequency dependent)."""

    def __cinit__(self, *args, **kwargs):
        self._rheology_ptr = NULL

    def __init__(self):
        cdef c_TideModelConfig config
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.Rheology, config)
        self._rheology_ptr = <c_RheologyTide*>ptr.get()
        self._tide_ptr     = move(ptr)
        self._ptr          = <c_TidalPyBaseClass*>self._tide_ptr.get()

    def __dealloc__(self):
        self._rheology_ptr = NULL


cdef class FixedQTide(TideBase):
    """Constant phase lag (fixed-Q): k_l(omega) = k_l * (1 - i / Q_l)."""

    def __cinit__(self, *args, **kwargs):
        self._fixedq_ptr = NULL

    def __init__(self, object fixed_k=None, object fixed_q=None):
        cdef c_TideModelConfig config
        config.fixed_k = _to_double_vector(fixed_k)
        config.fixed_q = _to_double_vector(fixed_q)
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.FixedQ, config)
        self._fixedq_ptr = <c_FixedQTide*>ptr.get()
        self._tide_ptr   = move(ptr)
        self._ptr        = <c_TidalPyBaseClass*>self._tide_ptr.get()

    def __dealloc__(self):
        self._fixedq_ptr = NULL

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree."""
        return self._fixedq_ptr.get_fixed_k(degree_l)

    def get_fixed_q(self, int degree_l) -> float:
        """Tidal quality factor Q_l at the given degree."""
        return self._fixedq_ptr.get_fixed_q(degree_l)

    cpdef dict get_config_dict(self):
        d = TideBase.get_config_dict(self)
        d["fixed_k"] = _per_degree_list(self.get_fixed_k)
        d["fixed_q"] = _per_degree_list(self.get_fixed_q)
        return d


cdef class FixedLagTide(TideBase):
    """Constant time lag (CTL): k_l(omega) = k_l * (1 - i * omega * dt_l)."""

    def __cinit__(self, *args, **kwargs):
        self._fixedlag_ptr = NULL

    def __init__(self, object fixed_k=None, object fixed_dt=None):
        cdef c_TideModelConfig config
        config.fixed_k  = _to_double_vector(fixed_k)
        config.fixed_dt = _to_double_vector(fixed_dt)
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.FixedLag, config)
        self._fixedlag_ptr = <c_FixedLagTide*>ptr.get()
        self._tide_ptr     = move(ptr)
        self._ptr          = <c_TidalPyBaseClass*>self._tide_ptr.get()

    def __dealloc__(self):
        self._fixedlag_ptr = NULL

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree."""
        return self._fixedlag_ptr.get_fixed_k(degree_l)

    def get_fixed_dt(self, int degree_l) -> float:
        """Tidal time lag dt_l [s] at the given degree."""
        return self._fixedlag_ptr.get_fixed_dt(degree_l)

    cpdef dict get_config_dict(self):
        d = TideBase.get_config_dict(self)
        d["fixed_k"]  = _per_degree_list(self.get_fixed_k)
        d["fixed_dt"] = _per_degree_list(self.get_fixed_dt)
        return d


cdef class CTLQTide(TideBase):
    """Constant time lag with a quality factor: k_l(omega) = k_l * (1 - i * omega * dt_l / Q_l)."""

    def __cinit__(self, *args, **kwargs):
        self._ctlq_ptr = NULL

    def __init__(self, object fixed_k=None, object fixed_dt=None, object fixed_q=None):
        cdef c_TideModelConfig config
        config.fixed_k  = _to_double_vector(fixed_k)
        config.fixed_dt = _to_double_vector(fixed_dt)
        config.fixed_q  = _to_double_vector(fixed_q)
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.CTLQ, config)
        self._ctlq_ptr = <c_CTLQTide*>ptr.get()
        self._tide_ptr = move(ptr)
        self._ptr      = <c_TidalPyBaseClass*>self._tide_ptr.get()

    def __dealloc__(self):
        self._ctlq_ptr = NULL

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree."""
        return self._ctlq_ptr.get_fixed_k(degree_l)

    def get_fixed_dt(self, int degree_l) -> float:
        """Tidal time lag dt_l [s] at the given degree."""
        return self._ctlq_ptr.get_fixed_dt(degree_l)

    def get_fixed_q(self, int degree_l) -> float:
        """Tidal quality factor Q_l at the given degree."""
        return self._ctlq_ptr.get_fixed_q(degree_l)

    cpdef dict get_config_dict(self):
        d = TideBase.get_config_dict(self)
        d["fixed_k"]  = _per_degree_list(self.get_fixed_k)
        d["fixed_dt"] = _per_degree_list(self.get_fixed_dt)
        d["fixed_q"]  = _per_degree_list(self.get_fixed_q)
        return d


# =====================================================================================================================
# Factory
# =====================================================================================================================
def make_tide(str model_name, dict config=None) -> TideBase:
    """Build a tide model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        One of 
          - ``"rheology"``
          - ``"cpl"`` / ``"fixed_q"``
          - ``"ctl"`` / ``"fixed_dt"``
          - ``"ctl_q"`` / ``"fixed_dt_q"``
    config : dict, optional
        Per-degree model parameters (``fixed_k``, ``fixed_q``, ``fixed_dt`` lists indexed
        from degree l = 2); absent keys default to zero / the C++ defaults.

    Returns
    -------
    TideBase
        The model subclass.

    Raises
    ------
    ValueError
        If the model name is unknown.
    """
    cdef c_TideModelConfig cfg = _build_tide_config(config)
    cdef c_TideModel model = c_tide_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_TideBase] ptr = c_find_tide(model, cfg)

    cdef RheologyTide rheology_tide
    cdef FixedQTide   fixedq_tide
    cdef FixedLagTide fixedlag_tide
    cdef CTLQTide     ctlq_tide
    if model == c_TideModel.Rheology:
        rheology_tide = RheologyTide.__new__(RheologyTide)
        rheology_tide._rheology_ptr = <c_RheologyTide*>ptr.get()
        rheology_tide._tide_ptr     = move(ptr)
        rheology_tide._ptr          = <c_TidalPyBaseClass*>rheology_tide._tide_ptr.get()
        return rheology_tide
    elif model == c_TideModel.FixedQ:
        fixedq_tide = FixedQTide.__new__(FixedQTide)
        fixedq_tide._fixedq_ptr = <c_FixedQTide*>ptr.get()
        fixedq_tide._tide_ptr   = move(ptr)
        fixedq_tide._ptr        = <c_TidalPyBaseClass*>fixedq_tide._tide_ptr.get()
        return fixedq_tide
    elif model == c_TideModel.FixedLag:
        fixedlag_tide = FixedLagTide.__new__(FixedLagTide)
        fixedlag_tide._fixedlag_ptr = <c_FixedLagTide*>ptr.get()
        fixedlag_tide._tide_ptr     = move(ptr)
        fixedlag_tide._ptr          = <c_TidalPyBaseClass*>fixedlag_tide._tide_ptr.get()
        return fixedlag_tide
    else:
        ctlq_tide = CTLQTide.__new__(CTLQTide)
        ctlq_tide._ctlq_ptr = <c_CTLQTide*>ptr.get()
        ctlq_tide._tide_ptr = move(ptr)
        ctlq_tide._ptr      = <c_TidalPyBaseClass*>ctlq_tide._tide_ptr.get()
        return ctlq_tide
