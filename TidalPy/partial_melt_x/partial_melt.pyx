# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's partial-melt models."""

import warnings

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass, cy_resolve_factory_config

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class PartialMeltBase(PhysicsBase):
    """Abstract base for partial-melt models. Instantiate a concrete subclass."""

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "PartialMeltBase is abstract; instantiate a concrete model "
            "(OffPartialMelt, SpohnPartialMelt, HenningPartialMelt).")

    def __dealloc__(self):
        self._melt_ptr.reset()
        self._ptr = NULL

    cdef void _adopt(self, unique_ptr[c_PartialMeltBase]& model) noexcept:
        """Take ownership of ``model``; the inherited ``_ptr`` observes it."""
        self._melt_ptr = move(model)
        self._ptr = <c_TidalPyBaseClass*>self._melt_ptr.get()

    @property
    def solidus(self) -> float:
        """Solidus temperature [K]."""
        self._check_ptr()
        return self._melt_ptr.get().get_solidus()

    @property
    def liquidus(self) -> float:
        """Liquidus temperature [K]."""
        self._check_ptr()
        return self._melt_ptr.get().get_liquidus()

    @property
    def liquid_shear(self) -> float:
        """Shear modulus of the fully-molten material, and the floor on the post-melt shear modulus [Pa]."""
        self._check_ptr()
        return self._melt_ptr.get().get_liquid_shear()

    @property
    def liquid_viscosity(self) -> float:
        """Viscosity of the fully-molten material, and the floor on the post-melt viscosity [Pa s]."""
        self._check_ptr()
        return self._melt_ptr.get().get_liquid_viscosity()

    @property
    def bulk_melt_weakening(self) -> bool:
        """Whether melt weakens the bulk modulus (``calc_bulk_modulus_melt``); off by default."""
        self._check_ptr()
        return self._melt_ptr.get().get_bulk_melt_weakening()

    @property
    def liquid_bulk_modulus(self) -> float:
        """Bulk modulus of the melt, used only when ``bulk_melt_weakening`` is on [Pa]."""
        self._check_ptr()
        return self._melt_ptr.get().get_liquid_bulk_modulus()

    def calc_melt_fraction(self, double temperature) -> float:
        """Volumetric melt fraction phi in [0, 1] from temperature [K]; NaN for a non-finite temperature."""
        self._check_ptr()
        return self._melt_ptr.get().calc_melt_fraction(temperature)

    def calc_partial_melt(
            self,
            double temperature,
            double premelt_viscosity,
            double premelt_shear) -> tuple:
        """Post-melt viscosity and shear modulus from the pre-melt state (all MKS).

        Both are floored at the model's liquid limits (``liquid_viscosity``, ``liquid_shear``). Below the solidus
        the pre-melt pair is returned; a non-finite temperature gives NaN.

        Returns
        -------
        (melt_fraction, postmelt_viscosity, postmelt_shear_modulus) : tuple of float
        """
        self._check_ptr()
        cdef c_PartialMeltInputs inputs
        inputs.temperature       = temperature
        inputs.premelt_viscosity = premelt_viscosity
        inputs.premelt_shear     = premelt_shear
        cdef c_PartialMeltResult result = self._melt_ptr.get().calc_partial_melt(inputs)
        return (result.melt_fraction, result.postmelt_viscosity, result.postmelt_shear_modulus)

    def calc_bulk_modulus_melt(
            self,
            double temperature,
            double premelt_bulk_modulus,
            double framework_shear_modulus) -> float:
        """Post-melt bulk modulus [Pa]; the pre-melt value unless ``bulk_melt_weakening`` is on.

        When on, the Hashin-Shtrikman (1963) bound for the melt (``liquid_bulk_modulus``) in a solid framework,
        evaluated with the framework's post-melt shear modulus ``framework_shear_modulus``: a weak reduction while
        the framework holds, the Reuss average of a suspension once its shear modulus has collapsed.
        """
        self._check_ptr()
        return self._melt_ptr.get().calc_bulk_modulus_melt(
            temperature, premelt_bulk_modulus, framework_shear_modulus)


cdef class OffPartialMelt(PartialMeltBase):
    """No melt weakening; post-melt strength equals pre-melt."""

    def __init__(
            self,
            double solidus=1600.0,
            double liquidus=2000.0,
            double liquid_shear=1.0e-5,
            double liquid_viscosity=0.2,
            cpp_bool bulk_melt_weakening=False,
            double liquid_bulk_modulus=2.0e10):
        cdef c_PartialMeltConfig config
        config.solidus             = solidus
        config.liquidus            = liquidus
        config.liquid_shear        = liquid_shear
        config.liquid_viscosity    = liquid_viscosity
        config.bulk_melt_weakening = bulk_melt_weakening
        config.liquid_bulk_modulus = liquid_bulk_modulus
        cdef unique_ptr[c_PartialMeltBase] model = c_find_partial_melt(c_PartialMeltModel.Off, config)
        self._adopt(model)


cdef class SpohnPartialMelt(PartialMeltBase):
    """Fischer & Spohn (1990) temperature-based melt law, anchored at the solidus.

    Above the solidus the post-melt viscosity and shear modulus are 10^(log10_at_solidus + s (1 / T - 1 / T_sol)),
    independent of the pre-melt values. The defaults reproduce Fischer and Spohn's fits, 10^(27000 / T - 1) Pa s and
    10^(82000 / T - 40.6) Pa, at the default 1600 K solidus.
    """

    def __init__(
            self,
            double solidus=1600.0,
            double liquidus=2000.0,
            double liquid_shear=1.0e-5,
            double fs_visc_power_slope=27000.0,
            double fs_visc_log10_at_solidus=15.875,
            double fs_shear_power_slope=82000.0,
            double fs_shear_log10_at_solidus=10.65,
            double liquid_viscosity=0.2,
            cpp_bool bulk_melt_weakening=False,
            double liquid_bulk_modulus=2.0e10):
        cdef c_PartialMeltConfig config
        config.solidus             = solidus
        config.liquidus            = liquidus
        config.liquid_shear        = liquid_shear
        config.liquid_viscosity    = liquid_viscosity
        config.bulk_melt_weakening = bulk_melt_weakening
        config.liquid_bulk_modulus = liquid_bulk_modulus
        config.fs_visc_power_slope       = fs_visc_power_slope
        config.fs_visc_log10_at_solidus  = fs_visc_log10_at_solidus
        config.fs_shear_power_slope      = fs_shear_power_slope
        config.fs_shear_log10_at_solidus = fs_shear_log10_at_solidus
        cdef unique_ptr[c_PartialMeltBase] model = c_find_partial_melt(c_PartialMeltModel.Spohn, config)
        self._adopt(model)

    @property
    def fs_visc_power_slope(self) -> float:
        """Viscosity-law temperature slope s [K] in 10^(log10_at_solidus + s (1 / T - 1 / T_sol)) [Pa s]."""
        self._check_ptr()
        return (<c_SpohnPartialMelt*>self._melt_ptr.get()).get_visc_power_slope()

    @property
    def fs_visc_log10_at_solidus(self) -> float:
        """log10 of the post-melt viscosity at the solidus [log10 Pa s]."""
        self._check_ptr()
        return (<c_SpohnPartialMelt*>self._melt_ptr.get()).get_visc_log10_at_solidus()

    @property
    def fs_shear_power_slope(self) -> float:
        """Shear-law temperature slope s [K] in 10^(log10_at_solidus + s (1 / T - 1 / T_sol)) [Pa]."""
        self._check_ptr()
        return (<c_SpohnPartialMelt*>self._melt_ptr.get()).get_shear_power_slope()

    @property
    def fs_shear_log10_at_solidus(self) -> float:
        """log10 of the post-melt shear modulus at the solidus [log10 Pa]."""
        self._check_ptr()
        return (<c_SpohnPartialMelt*>self._melt_ptr.get()).get_shear_log10_at_solidus()


cdef class HenningPartialMelt(PartialMeltBase):
    """Henning (2009/2010) three-regime melt weakening.

    Below the critical melt fraction the shear modulus is mu_pre exp[b1 (1 / T - 1 / T_sol)] with b1 =
    ``hn_shear_param_1``, which is 1 at the solidus. Henning et al. (2009) Eq. 20, exp(40000 / T - 25), is this law at
    the default 1600 K solidus.
    """

    def __init__(
            self,
            double solidus=1600.0,
            double liquidus=2000.0,
            double liquid_shear=1.0e-5,
            double crit_melt_frac=0.5,
            double crit_melt_frac_width=0.05,
            double hn_visc_slope_1=13.5,
            double hn_visc_falloff_slope=370.0,
            double hn_shear_param_1=40000.0,
            double hn_shear_falloff_slope=700.0,
            double liquid_viscosity=0.2,
            cpp_bool bulk_melt_weakening=False,
            double liquid_bulk_modulus=2.0e10):
        cdef c_PartialMeltConfig config
        config.solidus              = solidus
        config.liquidus             = liquidus
        config.liquid_shear         = liquid_shear
        config.liquid_viscosity     = liquid_viscosity
        config.bulk_melt_weakening  = bulk_melt_weakening
        config.liquid_bulk_modulus  = liquid_bulk_modulus
        config.crit_melt_frac       = crit_melt_frac
        config.crit_melt_frac_width = crit_melt_frac_width
        config.hn_visc_slope_1      = hn_visc_slope_1
        config.hn_visc_falloff_slope = hn_visc_falloff_slope
        config.hn_shear_param_1 = hn_shear_param_1
        config.hn_shear_falloff_slope = hn_shear_falloff_slope
        cdef unique_ptr[c_PartialMeltBase] model = c_find_partial_melt(c_PartialMeltModel.Henning, config)
        self._adopt(model)

    @property
    def crit_melt_frac(self) -> float:
        """Critical melt fraction phi_c at which the solid framework breaks down."""
        self._check_ptr()
        return (<c_HenningPartialMelt*>self._melt_ptr.get()).get_crit_melt_frac()

    @property
    def crit_melt_frac_width(self) -> float:
        """Width w of the transition band from phi_c to phi_c + w."""
        self._check_ptr()
        return (<c_HenningPartialMelt*>self._melt_ptr.get()).get_crit_melt_frac_width()

    @property
    def hn_visc_slope_1(self) -> float:
        """Viscosity weakening slope below phi_c: eta = eta_premelt * exp(-hn_visc_slope_1 * phi)."""
        self._check_ptr()
        return (<c_HenningPartialMelt*>self._melt_ptr.get()).get_visc_slope_1()

    @property
    def hn_visc_falloff_slope(self) -> float:
        """Viscosity falloff slope applied to (phi - phi_c) across the transition band."""
        self._check_ptr()
        return (<c_HenningPartialMelt*>self._melt_ptr.get()).get_visc_falloff_slope()

    @property
    def hn_shear_param_1(self) -> float:
        """Shear-law temperature parameter b_1 [K] in exp[b_1 (1 / T - 1 / T_sol)]."""
        self._check_ptr()
        return (<c_HenningPartialMelt*>self._melt_ptr.get()).get_shear_param_1()

    @property
    def hn_shear_falloff_slope(self) -> float:
        """Shear falloff slope applied to (phi - phi_c) across the transition band."""
        self._check_ptr()
        return (<c_HenningPartialMelt*>self._melt_ptr.get()).get_shear_falloff_slope()


# Every config key any partial-melt model reads; make_partial_melt rejects anything else.
PARTIAL_MELT_CONFIG_KEYS = frozenset({
    "solidus_k", "liquidus_k", "liquid_shear_pa", "liquid_viscosity_pas", "bulk_melt_weakening",
    "liquid_bulk_modulus_pa",
    "fs_visc_power_slope_k", "fs_visc_log10_at_solidus", "fs_shear_power_slope_k", "fs_shear_log10_at_solidus",
    "crit_melt_frac", "crit_melt_frac_width", "hn_visc_slope_1", "hn_visc_falloff_slope",
    "hn_shear_param_1_k", "hn_shear_falloff_slope"})

# Keys a 0.8.0 pre-release wrote into the user's TidalPy_Configs_x.toml that no model reads any more. They are dropped
# with a warning rather than rejected, so an existing configuration file still builds worlds.
RETIRED_PARTIAL_MELT_CONFIG_KEYS = {
    "hn_shear_param_2": "the Henning shear law is anchored at the solidus, exp[b1 (1/T - 1/T_sol)], so it has no "
                        "separate offset (the old default 25 is 40000 / 1600)",
}


# The wrapper class of each c_PartialMeltModel, in enum order.
_PARTIAL_MELT_CLASSES = (OffPartialMelt, SpohnPartialMelt, HenningPartialMelt)


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return (
        c_partial_melt_model_from_name(table_name.encode("utf-8"))
        == c_partial_melt_model_from_name(model_name.encode("utf-8")))


def make_partial_melt(str model_name, dict config=None) -> PartialMeltBase:
    """Build a partial-melt model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        ``"off"``/``"none"``, ``"spohn"``/``"fischer"``, or ``"henning"``.
    config : dict, optional
        Model parameters, keyed with their units (see ``PARTIAL_MELT_CONFIG_KEYS``). Absent keys fall back
        to the C++ defaults.

    Returns
    -------
    PartialMeltBase

    Raises
    ------
    ValueError
        Unknown model name, or a config key that no partial-melt model reads.

    Warns
    -----
    UserWarning
        A retired key (``RETIRED_PARTIAL_MELT_CONFIG_KEYS``) is present; it is ignored.
    """
    # The world builder's defaults hold no retired key (factory_defaults keeps only accepted ones), so only a
    # caller's config needs filtering.
    cdef list retired = [key for key in RETIRED_PARTIAL_MELT_CONFIG_KEYS if config and key in config]
    if retired:
        warnings.warn(
            "Partial-melt config key(s) no longer read and ignored: "
            + "; ".join(f"`{key}` ({RETIRED_PARTIAL_MELT_CONFIG_KEYS[key]})" for key in retired)
            + ". Remove them from your TidalPy_Configs_x.toml or world file.", stacklevel=2)
        config = {key: value for key, value in config.items() if key not in RETIRED_PARTIAL_MELT_CONFIG_KEYS}
    # None falls back to the same defaults the world-attached path uses.
    config = cy_resolve_factory_config(
        config, "material.partial_melt", PARTIAL_MELT_CONFIG_KEYS, model_name, _same_model, "partial-melt")
    # The default-constructed config carries the C++ defaults, so only override what the caller gave.
    cdef c_PartialMeltConfig cfg
    cfg.solidus                   = config.get("solidus_k", cfg.solidus)
    cfg.liquidus                  = config.get("liquidus_k", cfg.liquidus)
    cfg.liquid_shear              = config.get("liquid_shear_pa", cfg.liquid_shear)
    cfg.liquid_viscosity          = config.get("liquid_viscosity_pas", cfg.liquid_viscosity)
    cfg.bulk_melt_weakening       = bool(config.get("bulk_melt_weakening", cfg.bulk_melt_weakening))
    cfg.liquid_bulk_modulus       = config.get("liquid_bulk_modulus_pa", cfg.liquid_bulk_modulus)
    cfg.fs_visc_power_slope       = config.get("fs_visc_power_slope_k", cfg.fs_visc_power_slope)
    cfg.fs_visc_log10_at_solidus  = config.get("fs_visc_log10_at_solidus", cfg.fs_visc_log10_at_solidus)
    cfg.fs_shear_power_slope      = config.get("fs_shear_power_slope_k", cfg.fs_shear_power_slope)
    cfg.fs_shear_log10_at_solidus = config.get("fs_shear_log10_at_solidus", cfg.fs_shear_log10_at_solidus)
    cfg.crit_melt_frac            = config.get("crit_melt_frac", cfg.crit_melt_frac)
    cfg.crit_melt_frac_width      = config.get("crit_melt_frac_width", cfg.crit_melt_frac_width)
    cfg.hn_visc_slope_1           = config.get("hn_visc_slope_1", cfg.hn_visc_slope_1)
    cfg.hn_visc_falloff_slope     = config.get("hn_visc_falloff_slope", cfg.hn_visc_falloff_slope)
    cfg.hn_shear_param_1          = config.get("hn_shear_param_1_k", cfg.hn_shear_param_1)
    cfg.hn_shear_falloff_slope    = config.get("hn_shear_falloff_slope", cfg.hn_shear_falloff_slope)

    cdef c_PartialMeltModel model = c_partial_melt_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_PartialMeltBase] ptr = c_find_partial_melt(model, cfg)
    wrapper_class = _PARTIAL_MELT_CLASSES[<int>model]
    cdef PartialMeltBase wrapper = wrapper_class.__new__(wrapper_class)
    wrapper._adopt(ptr)
    return wrapper
