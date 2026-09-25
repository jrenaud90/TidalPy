# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's global (1D) tidal dissipation models.

A tide model turns a per-mode Love number into the dissipation multiplier -Im[k_l] the global mode collapse
uses. The per-degree fixed parameters are lists indexed from l = 2, and l = 2..10 is supported. The config
key carries the unit suffix; the constructor keyword and the C++ field stay ``fixed_dt``.
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
from TidalPy.Utilities_x.classes_x.classes import check_config_keys, factory_defaults
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Matches the C++ tide_.hpp slots.
cdef int C_TIDE_MIN_DEGREE = 2
cdef int C_TIDE_MAX_DEGREE = 10


cdef vector[double] cy_to_double_vector(object values) except *:
    """Convert an iterable of floats, indexed from l = 2, into a std::vector[double]."""
    cdef vector[double] out
    if values is None:
        return out
    for value in values:
        out.push_back(<double>value)
    return out


cdef c_TideModelConfig cy_build_tide_config(dict config) except *:
    """Build a c_TideModelConfig from a config dict with optional per-degree list keys."""
    cdef c_TideModelConfig cfg
    if config is None:
        return cfg
    if "fixed_k" in config:
        cfg.fixed_k = cy_to_double_vector(config["fixed_k"])
    if "fixed_q" in config:
        cfg.fixed_q = cy_to_double_vector(config["fixed_q"])
    if "fixed_dt_s" in config:
        cfg.fixed_dt = cy_to_double_vector(config["fixed_dt_s"])
    return cfg


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

    cdef void _adopt(self, unique_ptr[c_TideBase]& ptr) noexcept:
        """Take ownership of a built C++ model; `ptr` is left empty."""
        self._tide_ptr = move(ptr)
        self._ptr      = <c_TidalPyBaseClass*>self._tide_ptr.get()

    def calc_love_numbers(self, int degree_l, double frequency, LoveNumbers solver_love=None) -> LoveNumbers:
        """Full complex Love-number suite (k, h, l) at the tidal frequency [rad s-1].

        For analytic models h and l are NaN (no radial solution); the rheology model
        returns the supplied ``solver_love`` (a ``LoveNumbers``) unchanged. ``solver_love``
        is ignored by the analytic models and defaults to zeros when omitted.
        """
        self._check_ptr()
        cdef c_LoveNumbers solver_c
        if solver_love is not None:
            solver_c = solver_love._love
        cdef c_LoveNumbers result = self._tide_ptr.get().calc_love_numbers(degree_l, frequency, solver_c)
        cdef LoveNumbers out = LoveNumbers.__new__(LoveNumbers)
        out._love = result
        return out

    def calc_neg_imk(self, int degree_l, double frequency, LoveNumbers solver_love=None) -> float:
        """-Im[k_l] at the tidal frequency [rad s-1] (the mode-collapse dissipation multiplier)."""
        self._check_ptr()
        cdef c_LoveNumbers solver_c
        if solver_love is not None:
            solver_c = solver_love._love
        return self._tide_ptr.get().calc_neg_imk(degree_l, frequency, solver_c)

    @property
    def needs_radial_solve(self) -> bool:
        """Whether this model requires the radial solver to supply k_l (rheology)."""
        self._check_ptr()
        return bool(self._tide_ptr.get().needs_radial_solve())


cdef class RheologyTide(TideBase):
    """Rheology-based tide: k_l comes from the radial solver (frequency dependent)."""

    def __init__(self):
        cdef c_TideModelConfig config
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.Rheology, config)
        self._adopt(ptr)


cdef class FixedQTide(TideBase):
    """Constant phase lag (fixed-Q): k_l(omega) = k_l * (1 - i / Q_l)."""

    def __init__(self, object fixed_k=None, object fixed_q=None):
        cdef c_TideModelConfig config
        config.fixed_k = cy_to_double_vector(fixed_k)
        config.fixed_q = cy_to_double_vector(fixed_q)
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.FixedQ, config)
        self._adopt(ptr)

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree."""
        self._check_ptr()
        return (<c_FixedQTide*>self._tide_ptr.get()).get_fixed_k(degree_l)

    def get_fixed_q(self, int degree_l) -> float:
        """Tidal quality factor Q_l at the given degree."""
        self._check_ptr()
        return (<c_FixedQTide*>self._tide_ptr.get()).get_fixed_q(degree_l)


cdef class FixedLagTide(TideBase):
    """Constant time lag (CTL): k_l(omega) = k_l * (1 - i * omega * dt_l)."""

    def __init__(self, object fixed_k=None, object fixed_dt=None):
        cdef c_TideModelConfig config
        config.fixed_k  = cy_to_double_vector(fixed_k)
        config.fixed_dt = cy_to_double_vector(fixed_dt)
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.FixedLag, config)
        self._adopt(ptr)

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree."""
        self._check_ptr()
        return (<c_FixedLagTide*>self._tide_ptr.get()).get_fixed_k(degree_l)

    def get_fixed_dt(self, int degree_l) -> float:
        """Tidal time lag dt_l [s] at the given degree."""
        self._check_ptr()
        return (<c_FixedLagTide*>self._tide_ptr.get()).get_fixed_dt(degree_l)


cdef class CTLQTide(TideBase):
    """Constant time lag with a quality factor: k_l(omega) = k_l * (1 - i * omega * dt_l / Q_l)."""

    def __init__(self, object fixed_k=None, object fixed_dt=None, object fixed_q=None):
        cdef c_TideModelConfig config
        config.fixed_k  = cy_to_double_vector(fixed_k)
        config.fixed_dt = cy_to_double_vector(fixed_dt)
        config.fixed_q  = cy_to_double_vector(fixed_q)
        cdef unique_ptr[c_TideBase] ptr = c_find_tide(c_TideModel.CTLQ, config)
        self._adopt(ptr)

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree."""
        self._check_ptr()
        return (<c_CTLQTide*>self._tide_ptr.get()).get_fixed_k(degree_l)

    def get_fixed_dt(self, int degree_l) -> float:
        """Tidal time lag dt_l [s] at the given degree."""
        self._check_ptr()
        return (<c_CTLQTide*>self._tide_ptr.get()).get_fixed_dt(degree_l)

    def get_fixed_q(self, int degree_l) -> float:
        """Tidal quality factor Q_l at the given degree."""
        self._check_ptr()
        return (<c_CTLQTide*>self._tide_ptr.get()).get_fixed_q(degree_l)


# The wrapper class of each c_TideModel, indexed by its enum value.
cdef tuple cy_tide_classes = (RheologyTide, FixedQTide, FixedLagTide, CTLQTide)


# Every config key any tide model reads; make_tide rejects anything else.
TIDE_CONFIG_KEYS = frozenset({"fixed_k", "fixed_q", "fixed_dt_s"})


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return c_tide_model_from_name(table_name.lower().encode("utf-8")) == c_tide_model_from_name(model_name.lower().encode("utf-8"))


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
        Per-degree parameters (``fixed_k``, ``fixed_q``, ``fixed_dt_s`` [s]), indexed from l = 2. A key left out
        takes the ``[tides]`` value of the TidalPy configuration, as the world builder does, so
        ``{"fixed_q": [50]}`` keeps the configured ``fixed_k``.

    Returns
    -------
    TideBase

    Raises
    ------
    ValueError
        Unknown model name, or a config key outside those three.
    """
    # Key by key over the same defaults the world-attached path uses, so a partial config never leaves a list empty
    # (an empty fixed_k would silently give no dissipation).
    if config is not None:
        check_config_keys(config, TIDE_CONFIG_KEYS, "tide")
    config = {**factory_defaults("tides", TIDE_CONFIG_KEYS, model_name, _same_model), **(config or {})}
    cdef c_TideModelConfig cfg = cy_build_tide_config(config)
    cdef c_TideModel model = c_tide_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_TideBase] ptr = c_find_tide(model, cfg)
    tide_class = cy_tide_classes[<int>model]
    cdef TideBase tide = tide_class.__new__(tide_class)
    tide._adopt(ptr)
    return tide
