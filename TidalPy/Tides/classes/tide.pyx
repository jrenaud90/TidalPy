# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's global (1D) tidal dissipation models.

A tide model turns a per-mode Love number into the dissipation multiplier -Im[k_l] the global mode collapse
uses. Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table. The analytic
models' parameters are lists indexed from l = 2, and l = 2..10 is supported; a list left out takes its ``[tides]``
value of the TidalPy configuration.
"""

from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_share_physics,
    cy_collect_parameters,
    cy_param_map,
)
from TidalPy.Utilities.classes.classes import factory_defaults
from TidalPy.Utilities.classes.families import ModelFamily
from TidalPy.Tides.love.love cimport LoveNumbers, c_LoveNumbers

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class TideBase(PhysicsBase):
    """Base for global tide models, which give the per-mode Love number of the global mode collapse.

    Instantiate a concrete model (``RheologyTide``, ``FixedQTide``, ``FixedLagTide``, ``CTLQTide``) with its
    parameters positionally (in the order ``get_parameter_info()`` lists them) or as keywords (argument names or config
    keys), or build one by name with ``make_tide``. A per-degree list left out (or None) takes its ``[tides]`` value of
    the TidalPy configuration, so ``FixedQTide([0.3])`` has the configured Q rather than none. A world shares the model
    it is given (``set_tide_model``).
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "TideBase is abstract; instantiate a concrete model (RheologyTide, FixedQTide, FixedLagTide, "
                "CTLQTide) or call make_tide.")
        cdef dict merged = cy_collect_parameters(type(self), args, config, parameters)
        for key, value in factory_defaults("tides", _FAMILY.config_keys_of(model_name)).items():
            if merged.get(key) is None:
                merged[key] = value
        cdef unique_ptr[c_TideBase] model = c_find_tide((<str>model_name).encode("utf-8"), cy_param_map(merged))
        self._set_model(c_share_physics[c_TideBase](move(model)))

    @classmethod
    def _default_model(cls):
        """The model at its own defaults (empty lists), which reads no configuration (for its parameter
        descriptions)."""
        cdef TideBase model = cls.__new__(cls)
        cdef c_ParamMap no_parameters
        cdef unique_ptr[c_TideBase] built = c_find_tide((<str>cls.MODEL_NAME).encode("utf-8"), no_parameters)
        model._set_model(c_share_physics[c_TideBase](move(built)))
        return model

    cdef c_TideBase* _tide(self) except NULL:
        self._check_ptr()
        return <c_TideBase*>self._ptr

    def calc_love_numbers(self, int degree_l, double frequency, LoveNumbers solver_love=None) -> LoveNumbers:
        """Full complex Love-number suite (k, h, l) at the tidal frequency [rad s-1].

        For analytic models h and l are NaN (no radial solution); the rheology model
        returns the supplied ``solver_love`` (a ``LoveNumbers``) unchanged. ``solver_love``
        is ignored by the analytic models and defaults to zeros when omitted.
        """
        cdef c_LoveNumbers solver_c
        if solver_love is not None:
            solver_c = solver_love._love
        cdef c_LoveNumbers result = self._tide().calc_love_numbers(degree_l, frequency, solver_c)
        cdef LoveNumbers out = LoveNumbers.__new__(LoveNumbers)
        out._love = result
        return out

    def calc_neg_imk(self, int degree_l, double frequency, LoveNumbers solver_love=None) -> float:
        """-Im[k_l] at the tidal frequency [rad s-1] (the mode-collapse dissipation multiplier)."""
        cdef c_LoveNumbers solver_c
        if solver_love is not None:
            solver_c = solver_love._love
        return self._tide().calc_neg_imk(degree_l, frequency, solver_c)

    def get_fixed_k(self, int degree_l) -> float:
        """Static potential Love number k_l at the given degree (0 past the end of the list); NaN when the model
        carries none."""
        return self._tide().get_fixed_k(degree_l)

    def get_fixed_q(self, int degree_l) -> float:
        """Tidal quality factor Q_l at the given degree (0 past the end of the list); NaN when the model carries
        none."""
        return self._tide().get_fixed_q(degree_l)

    def get_fixed_dt(self, int degree_l) -> float:
        """Tidal time lag dt_l [s] at the given degree (0 past the end of the list); NaN when the model carries
        none."""
        return self._tide().get_fixed_dt(degree_l)

    @property
    def needs_radial_solve(self) -> bool:
        """Whether this model requires the radial solver to supply k_l (rheology)."""
        return True if self._tide().needs_radial_solve() else False


cdef class RheologyTide(TideBase):
    """Rheology-based tide: k_l comes from the radial solver (frequency dependent). Takes no parameters."""
    MODEL_NAME = "rheology"


cdef class FixedQTide(TideBase):
    """Constant phase lag (aliases ``"cpl"``, ``"constant_phase_lag"``): k_l(omega) = k_l * (1 - i / Q_l)."""
    MODEL_NAME = "fixed_q"


cdef class FixedLagTide(TideBase):
    """Constant time lag (aliases ``"ctl"``, ``"constant_time_lag"``): k_l(omega) = k_l * (1 - i * |omega| * dt_l)."""
    MODEL_NAME = "fixed_dt"


cdef class CTLQTide(TideBase):
    """Constant time lag with a quality factor (aliases ``"ctl_q"``, ``"constant_time_lag_and_q"``):
    k_l(omega) = k_l * (1 - i * |omega| * dt_l / Q_l)."""
    MODEL_NAME = "fixed_dt_q"


# The family's name lookup: the C++ registry's alias-aware, case-insensitive canonical name.
_FAMILY = ModelFamily(
    "tide",
    (RheologyTide, FixedQTide, FixedLagTide, CTLQTide),
    lambda model_name: c_tide_canonical_name(model_name.encode("utf-8")).decode("utf-8"))

# Every config key any tide model reads.
TIDE_CONFIG_KEYS = _FAMILY.config_keys


def tide_model_names() -> tuple:
    """The canonical names of the tide models."""
    return _FAMILY.model_names()


def tide_config_keys(str model_name) -> frozenset:
    """The config keys a tide model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _FAMILY.config_keys_of(model_name)


def make_tide(str model_name, dict config=None, **parameters) -> TideBase:
    """Build a tide model from a (case-insensitive) name, a config dict, and keyword parameters.

    Parameters
    ----------
    model_name : str
        ``rheology``, ``fixed_q`` (``cpl``), ``fixed_dt`` (``ctl``), or ``fixed_dt_q`` (``ctl_q``).
    config : dict, optional
        Model parameters by config key (``fixed_k``, ``fixed_q``, ``fixed_dt_s`` [s], lists indexed from l = 2, as
        the model reads them; see ``get_parameter_info()``). A list left out takes the ``[tides]`` value of the
        TidalPy configuration, as the world builder does, so ``{"fixed_q": [50]}`` keeps the configured ``fixed_k``.
    **parameters
        The same parameters by argument name or config key (``fixed_q=[50]``, ``fixed_dt=[600.0]``), merged over
        ``config``.

    Returns
    -------
    TideBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    TypeError
        One parameter given under both of its spellings.
    """
    return _FAMILY.make(model_name, config, **parameters)
