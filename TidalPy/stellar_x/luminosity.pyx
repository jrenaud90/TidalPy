# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's stellar luminosity models.

References
----------
- Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 - low-mass mass-luminosity polynomial exponent.
"""

from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp

# The shared vector helpers build their result arrays through the NumPy C API.
cnp.import_array()

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass
from TidalPy.Utilities_x.classes_x.classes import check_config_keys

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef object cy_solve_luminosity(c_LuminosityBase* model, object mass):
    """Luminosity for a float or ndarray mass; a float for a float mass, else an array of the mass's shape."""
    cdef vector[vector[double]] inputs
    cdef vector[double] luminosity
    cdef object shape = cy_broadcast_inputs((mass,), inputs, False)
    if shape is None:
        return float(model.calc_luminosity(<double>mass))
    with nogil:
        model.calc_luminosity_vectorize_mass(inputs[0], luminosity)
    return cy_vector_to_ndarray(luminosity, shape)


cdef class LuminosityBase(PhysicsBase):
    """Abstract base for stellar luminosity models; owns the most-derived C++ model object."""

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "LuminosityBase is abstract; instantiate a concrete model "
            "(FixedLuminosity, MassToLuminosity, PowerLawLuminosity)."
        )

    def __dealloc__(self):
        # unique_ptr frees the most-derived C++ object; _ptr is only an observer.
        self._luminosity_ptr.reset()
        self._ptr = NULL

    cdef void _adopt(self, unique_ptr[c_LuminosityBase]& model) noexcept:
        """Take ownership of ``model``; the inherited ``_ptr`` observes it."""
        self._luminosity_ptr = move(model)
        self._ptr = <c_TidalPyBaseClass*>self._luminosity_ptr.get()

    def calc_luminosity(self, mass):
        """Stellar luminosity [W] from mass.

        Parameters
        ----------
        mass : float or numpy.ndarray
            Stellar mass [kg].

        Returns
        -------
        float or numpy.ndarray
            Stellar luminosity [W]; a float for scalar mass, else a same-shape float64 array.

        Notes
        -----
        Assumes main-sequence mass-luminosity scaling.
        """
        self._check_ptr()
        return cy_solve_luminosity(self._luminosity_ptr.get(), mass)

    def calc_luminosity_from_temperature(self, double temperature, double radius) -> float:
        """Stefan-Boltzmann luminosity [W] = 4*pi*R^2*sigma*T^4; NaN for non-positive inputs."""
        self._check_ptr()
        return self._luminosity_ptr.get().calc_luminosity_from_temperature(temperature, radius)

    def calc_temperature_from_luminosity(self, double luminosity, double radius) -> float:
        """Effective temperature ``T = (L / (4*pi*R^2*sigma))^(1/4)`` [K]; NaN for non-positive inputs."""
        self._check_ptr()
        return self._luminosity_ptr.get().calc_temperature_from_luminosity(luminosity, radius)

    def calc_effective_temperature(self, double mass, double radius) -> float:
        """Effective temperature [K] derived from the stellar mass (mass -> L -> T)."""
        self._check_ptr()
        return self._luminosity_ptr.get().calc_effective_temperature(mass, radius)


cdef class FixedLuminosity(LuminosityBase):
    """Luminosity supplied directly (mass independent).

    Parameters
    ----------
    luminosity : float, optional
        The luminosity [W] to report regardless of mass. Default ``0.0``.
    """

    def __init__(self, double luminosity=0.0):
        cdef c_LuminosityConfig config
        config.luminosity = luminosity
        cdef unique_ptr[c_LuminosityBase] model = c_find_luminosity(c_LuminosityModel.Fixed, config)
        self._adopt(model)

    @property
    def luminosity(self) -> float:
        """The stored luminosity [W]."""
        self._check_ptr()
        return (<c_FixedLuminosity*>self._luminosity_ptr.get()).get_luminosity()


cdef class MassToLuminosity(LuminosityBase):
    """Piecewise main-sequence L(M): the Cuntz and Wang (2018) polynomial exponent at low mass, the
    standard power-law regimes elsewhere. Takes no parameters.
    """

    def __init__(self):
        cdef c_LuminosityConfig config
        cdef unique_ptr[c_LuminosityBase] model = c_find_luminosity(c_LuminosityModel.MassToLuminosity, config)
        self._adopt(model)


cdef class PowerLawLuminosity(LuminosityBase):
    """Single power law ``L = Lsun * coeff * (M / Msun)^exponent``.

    Parameters
    ----------
    coeff : float, optional
        Dimensionless prefactor. Default ``1.0``.
    exponent : float, optional
        Dimensionless exponent. Default ``3.5`` (classic main-sequence value).
    """

    def __init__(self, double coeff=1.0, double exponent=3.5):
        cdef c_LuminosityConfig config
        config.power_law_coeff    = coeff
        config.power_law_exponent = exponent
        cdef unique_ptr[c_LuminosityBase] model = c_find_luminosity(c_LuminosityModel.PowerLaw, config)
        self._adopt(model)

    @property
    def coeff(self) -> float:
        """Dimensionless prefactor."""
        self._check_ptr()
        return (<c_PowerLawLuminosity*>self._luminosity_ptr.get()).get_coeff()

    @property
    def exponent(self) -> float:
        """Dimensionless exponent."""
        self._check_ptr()
        return (<c_PowerLawLuminosity*>self._luminosity_ptr.get()).get_exponent()


# Every config key any luminosity model reads; make_luminosity rejects anything else.
LUMINOSITY_CONFIG_KEYS = frozenset({"luminosity_w", "power_law_coeff", "power_law_exponent"})


# The wrapper class of each c_LuminosityModel, in enum order.
_LUMINOSITY_CLASSES = (FixedLuminosity, MassToLuminosity, PowerLawLuminosity)


def make_luminosity(str model_name, dict config=None):
    """Build a luminosity model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        ``fixed`` (``constant``), ``mass_to_luminosity`` (``cuntz_wang``, ``cw``), or ``power_law``.
    config : dict, optional
        Model parameters. For ``fixed``: ``luminosity_w``. For ``power_law``: ``power_law_coeff``,
        ``power_law_exponent``. ``mass_to_luminosity`` takes none.

    Returns
    -------
    LuminosityBase

    Raises
    ------
    ValueError
        Unknown model name, or a config key that no luminosity model reads.
    """
    check_config_keys(config, LUMINOSITY_CONFIG_KEYS, "luminosity")
    if config is None:
        config = {}

    # The default-constructed config carries the C++ defaults, so only override what the caller gave.
    cdef c_LuminosityConfig cfg
    cfg.luminosity         = config.get("luminosity_w", cfg.luminosity)
    cfg.power_law_coeff    = config.get("power_law_coeff", cfg.power_law_coeff)
    cfg.power_law_exponent = config.get("power_law_exponent", cfg.power_law_exponent)

    cdef c_LuminosityModel model = c_luminosity_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_LuminosityBase] ptr = c_find_luminosity(model, cfg)
    wrapper_class = _LUMINOSITY_CLASSES[<int>model]
    cdef LuminosityBase wrapper = wrapper_class.__new__(wrapper_class)
    wrapper._adopt(ptr)
    return wrapper


# Convenience functions. Each builds a stack-allocated C++ model that dies with the call.

def fixed(mass, double luminosity=0.0):
    """Luminosity for the Fixed model [W] (returns ``luminosity`` regardless of mass)."""
    cdef c_LuminosityConfig cfg
    cfg.luminosity = luminosity
    cdef c_FixedLuminosity model = c_FixedLuminosity(cfg)
    return cy_solve_luminosity(<c_LuminosityBase*>&model, mass)


def mass_to_luminosity(mass):
    """Luminosity for the MassToLuminosity model [W] (piecewise main-sequence relation)."""
    cdef c_LuminosityConfig cfg
    cdef c_MassToLuminosity model = c_MassToLuminosity(cfg)
    return cy_solve_luminosity(<c_LuminosityBase*>&model, mass)


def power_law(mass, double coeff=1.0, double exponent=3.5):
    """Luminosity for the PowerLaw model [W]: ``L = Lsun * coeff * (M / Msun)^exponent``."""
    cdef c_LuminosityConfig cfg
    cfg.power_law_coeff    = coeff
    cfg.power_law_exponent = exponent
    cdef c_PowerLawLuminosity model = c_PowerLawLuminosity(cfg)
    return cy_solve_luminosity(<c_LuminosityBase*>&model, mass)
