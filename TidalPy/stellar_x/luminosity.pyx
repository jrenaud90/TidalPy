# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's stellar luminosity models.

- ``FixedLuminosity``    (alias ``"constant"``)              - luminosity set directly (mass independent).
- ``MassToLuminosity``   (aliases ``"cuntz_wang"``/``"cw"``) - piecewise main-sequence L(M) relation.
- ``PowerLawLuminosity`` (alias ``"power_law"``)             - single power law L = Lsun*coeff*(M/Msun)^p.

Each model returns a star's luminosity [W] from its mass [kg] and shares the Stefan-Boltzmann
effective-temperature conversions.

References
----------
- Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 - low-mass mass-luminosity polynomial exponent.
"""

from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector

import numpy as np

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass
from TidalPy.Utilities_x.classes_x.classes import check_config_keys

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# Internal helpers for vectorized solving
# =====================================================================================================================
cdef void _fill_vector(double[::1] src, vector[double]& dst) noexcept:
    """Copy a contiguous 1-D float64 memoryview into a std::vector[double]."""
    cdef Py_ssize_t n = src.shape[0]
    cdef Py_ssize_t i
    dst.resize(n)
    for i in range(n):
        dst[i] = src[i]


cdef object _double_vector_to_ndarray(vector[double]& src, tuple shape):
    """Build a float64 ndarray (of the given shape) from a std::vector."""
    cdef Py_ssize_t n = <Py_ssize_t>src.size()
    cdef Py_ssize_t i
    out = np.empty(n, dtype=np.float64)
    cdef double[::1] mv = out
    for i in range(n):
        mv[i] = src[i]
    return out.reshape(shape)


cdef object _solve_luminosity(c_LuminosityBase* model, object mass):
    """Solve luminosity for float and/or np.ndarray mass inputs.

    Returns a Python ``float`` for scalar mass, else a float64 ``np.ndarray`` (same shape).
    """
    cdef vector[double] vmass
    cdef vector[double] vout
    cdef double[::1] mv

    if not isinstance(mass, np.ndarray):
        return float(model.calc_luminosity(<double>mass))

    mass_arr = np.ascontiguousarray(mass, dtype=np.float64)
    mv = mass_arr.ravel()
    _fill_vector(mv, vmass)
    model.calc_luminosity_vectorize_mass(vmass, vout)
    return _double_vector_to_ndarray(vout, mass_arr.shape)


# =====================================================================================================================
# LuminosityBase
# =====================================================================================================================
cdef class LuminosityBase(PhysicsBase):
    """Abstract base for all stellar luminosity models.

    Holds the owning ``unique_ptr`` to the most-derived C++ model object and exposes the shared
    luminosity calculations (mass -> luminosity and the Stefan-Boltzmann temperature conversions). Not
    directly instantiable - construct a concrete model (e.g. ``MassToLuminosity``).
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_LuminosityBase> auto-inits to nullptr; concrete models set it

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "LuminosityBase is abstract; instantiate a concrete model "
            "(FixedLuminosity, MassToLuminosity, PowerLawLuminosity)."
        )

    def __dealloc__(self):
        # unique_ptr frees the most-derived C++ object; null the base observer ptr.
        self._luminosity_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_luminosity(self, mass):
        """Stellar luminosity [W] from mass.

        Parameters
        ----------
        mass : float or numpy.ndarray
            Stellar mass [kg].

        Returns
        -------
        float or numpy.ndarray
            Stellar luminosity [W] (float for scalar mass, else same-shape float64 array).

        Assumptions
        -----------
        - Main-sequence mass-luminosity scaling (model specific).
        """
        self._check_ptr()
        return _solve_luminosity(self._luminosity_ptr.get(), mass)

    def calc_luminosity_from_temperature(self, double temperature, double radius) -> float:
        """Stefan-Boltzmann luminosity [W] = 4*pi*R^2*sigma*T^4.

        Returns NaN for a non-positive temperature/radius.
        """
        self._check_ptr()
        return self._luminosity_ptr.get().calc_luminosity_from_temperature(temperature, radius)

    def calc_temperature_from_luminosity(self, double luminosity, double radius) -> float:
        """Effective temperature [K] from luminosity via Stefan-Boltzmann.

        ``T = (L / (4*pi*R^2*sigma))^(1/4)``. Returns NaN for a non-positive luminosity/radius.
        """
        self._check_ptr()
        return self._luminosity_ptr.get().calc_temperature_from_luminosity(luminosity, radius)

    def calc_effective_temperature(self, double mass, double radius) -> float:
        """Effective temperature [K] derived from the stellar mass (mass -> L -> T)."""
        self._check_ptr()
        return self._luminosity_ptr.get().calc_effective_temperature(mass, radius)


# =====================================================================================================================
# FixedLuminosity (alias "constant")
# =====================================================================================================================
cdef class FixedLuminosity(LuminosityBase):
    """Luminosity supplied directly (mass independent).

    Parameters
    ----------
    luminosity : float, optional
        The luminosity [W] to report regardless of mass. Default ``0.0``.
    """

    def __cinit__(self, *args, **kwargs):
        self._fixed_ptr = NULL

    def __init__(self, double luminosity=0.0):
        cdef c_LuminosityConfig config
        config.luminosity = luminosity
        # Build through the C++ factory (make_unique) and adopt ownership; no raw new/delete.
        cdef unique_ptr[c_LuminosityBase] ptr = c_find_luminosity(c_LuminosityModel.Fixed, config)
        self._fixed_ptr = <c_FixedLuminosity*>ptr.get()
        self._luminosity_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._luminosity_ptr.get()

    def __dealloc__(self):
        self._fixed_ptr = NULL  # base unique_ptr owns the C++ object

    @property
    def luminosity(self) -> float:
        """The stored luminosity [W]."""
        self._check_ptr()
        return self._fixed_ptr.get_luminosity()


# =====================================================================================================================
# MassToLuminosity (aliases "cuntz_wang", "cw")
# =====================================================================================================================
cdef class MassToLuminosity(LuminosityBase):
    """Piecewise main-sequence mass-luminosity relation.

    Uses the Cuntz and Wang (2018) polynomial exponent for low-mass stars and the standard
    main-sequence power-law regimes elsewhere. Takes no configurable parameters.
    """

    def __init__(self):
        cdef c_LuminosityConfig config
        cdef unique_ptr[c_LuminosityBase] ptr = c_find_luminosity(c_LuminosityModel.MassToLuminosity, config)
        self._luminosity_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._luminosity_ptr.get()


# =====================================================================================================================
# PowerLawLuminosity (alias "power_law")
# =====================================================================================================================
cdef class PowerLawLuminosity(LuminosityBase):
    """Single power-law mass-luminosity relation.

    ``L = Lsun * coeff * (M / Msun)^exponent``.

    Parameters
    ----------
    coeff : float, optional
        Dimensionless prefactor. Default ``1.0``.
    exponent : float, optional
        Dimensionless exponent. Default ``3.5`` (classic main-sequence value).
    """

    def __cinit__(self, *args, **kwargs):
        self._power_law_ptr = NULL

    def __init__(self, double coeff=1.0, double exponent=3.5):
        cdef c_LuminosityConfig config
        config.power_law_coeff    = coeff
        config.power_law_exponent = exponent
        cdef unique_ptr[c_LuminosityBase] ptr = c_find_luminosity(c_LuminosityModel.PowerLaw, config)
        self._power_law_ptr = <c_PowerLawLuminosity*>ptr.get()
        self._luminosity_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._luminosity_ptr.get()

    def __dealloc__(self):
        self._power_law_ptr = NULL

    @property
    def coeff(self) -> float:
        """Dimensionless prefactor."""
        self._check_ptr()
        return self._power_law_ptr.get_coeff()

    @property
    def exponent(self) -> float:
        """Dimensionless exponent."""
        self._check_ptr()
        return self._power_law_ptr.get_exponent()


# =====================================================================================================================
# Factory
# =====================================================================================================================
# Every config key some luminosity model reads; make_luminosity rejects anything else.
LUMINOSITY_CONFIG_KEYS = frozenset({"luminosity_w", "power_law_coeff", "power_law_exponent"})


def make_luminosity(str model_name, dict config=None):
    """Build a luminosity model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias. Recognized names: ``fixed`` (``constant``), ``mass_to_luminosity``
        (``cuntz_wang``, ``cw``), ``power_law`` (``powerlaw``).
    config : dict, optional
        Model parameters. For ``fixed``: ``luminosity_w``. For ``power_law``: ``power_law_coeff``,
        ``power_law_exponent``. ``mass_to_luminosity`` takes no parameters.

    Returns
    -------
    LuminosityBase
        A concrete luminosity model instance.

    Raises
    ------
    ValueError
        If the model name is not recognized, or if ``config`` holds a key that no luminosity model reads.
    """
    check_config_keys(config, LUMINOSITY_CONFIG_KEYS, "luminosity")
    if config is None:
        config = {}

    # A default-constructed config carries the C++ defaults; override only what the caller supplies.
    cdef c_LuminosityConfig cfg
    if "luminosity_w" in config:
        cfg.luminosity = config["luminosity_w"]
    if "power_law_coeff" in config:
        cfg.power_law_coeff = config["power_law_coeff"]
    if "power_law_exponent" in config:
        cfg.power_law_exponent = config["power_law_exponent"]

    # Map name/alias to the enum (ValueError on an unknown name through except +).
    cdef c_LuminosityModel model = c_luminosity_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_LuminosityBase] ptr = c_find_luminosity(model, cfg)

    # Adopt the owning unique_ptr into the matching rich Python wrapper.
    cdef FixedLuminosity    fixed_wrapper
    cdef MassToLuminosity   mass_wrapper
    cdef PowerLawLuminosity power_wrapper

    if model == c_LuminosityModel.Fixed:
        fixed_wrapper = FixedLuminosity.__new__(FixedLuminosity)
        fixed_wrapper._fixed_ptr      = <c_FixedLuminosity*>ptr.get()
        fixed_wrapper._luminosity_ptr = move(ptr)
        fixed_wrapper._ptr            = <c_TidalPyBaseClass*>fixed_wrapper._fixed_ptr
        return fixed_wrapper
    elif model == c_LuminosityModel.MassToLuminosity:
        mass_wrapper = MassToLuminosity.__new__(MassToLuminosity)
        mass_wrapper._luminosity_ptr = move(ptr)
        mass_wrapper._ptr            = <c_TidalPyBaseClass*>mass_wrapper._luminosity_ptr.get()
        return mass_wrapper
    else:  # c_LuminosityModel.PowerLaw
        power_wrapper = PowerLawLuminosity.__new__(PowerLawLuminosity)
        power_wrapper._power_law_ptr  = <c_PowerLawLuminosity*>ptr.get()
        power_wrapper._luminosity_ptr = move(ptr)
        power_wrapper._ptr            = <c_TidalPyBaseClass*>power_wrapper._power_law_ptr
        return power_wrapper


# =====================================================================================================================
# Direct luminosity convenience functions
#
# Each builds a stack-allocated C++ model, solves, and returns a float for scalar mass or a
# same-shape float64 ndarray.
# =====================================================================================================================

def fixed(mass, double luminosity=0.0):
    """Luminosity for the Fixed model [W] (returns ``luminosity`` regardless of mass)."""
    cdef c_LuminosityConfig cfg
    cfg.luminosity = luminosity
    cdef c_FixedLuminosity model = c_FixedLuminosity(cfg)
    return _solve_luminosity(<c_LuminosityBase*>&model, mass)


def mass_to_luminosity(mass):
    """Luminosity for the MassToLuminosity model [W] (piecewise main-sequence relation)."""
    cdef c_LuminosityConfig cfg
    cdef c_MassToLuminosity model = c_MassToLuminosity(cfg)
    return _solve_luminosity(<c_LuminosityBase*>&model, mass)


def power_law(mass, double coeff=1.0, double exponent=3.5):
    """Luminosity for the PowerLaw model [W]: ``L = Lsun * coeff * (M / Msun)^exponent``."""
    cdef c_LuminosityConfig cfg
    cfg.power_law_coeff    = coeff
    cfg.power_law_exponent = exponent
    cdef c_PowerLawLuminosity model = c_PowerLawLuminosity(cfg)
    return _solve_luminosity(<c_LuminosityBase*>&model, mass)
