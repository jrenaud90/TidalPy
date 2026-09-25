# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's rheology models. Each returns a complex modulus mu* [Pa].

References
----------
- Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000
- Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
- Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784
"""

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.memory cimport unique_ptr
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
from TidalPy.Utilities_x.arrays.vectors cimport cy_broadcast_inputs, cy_complex_vector_to_ndarray
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass, cy_resolve_factory_config

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef object cy_solve_complex_modulus(
        c_RheologyBase* model, object modulus, object viscosity, object frequency, cpp_bool flatten):
    """Complex modulus for float or ndarray inputs broadcast together (cy_broadcast_inputs); a Python complex when
    every input is a float and ``flatten`` is off."""
    cdef cpp_complex[double] scalar_result
    cdef vector[vector[double]] inputs
    cdef vector[cpp_complex[double]] complex_modulus
    cdef object shape = cy_broadcast_inputs((modulus, viscosity, frequency), inputs, flatten)
    if shape is None:
        scalar_result = model.calc_complex_modulus(<double>modulus, <double>viscosity, <double>frequency)
        return complex(scalar_result.real(), scalar_result.imag())
    with nogil:
        model.calc_complex_modulus_vectorize(inputs[0], inputs[1], inputs[2], complex_modulus)
    return cy_complex_vector_to_ndarray(complex_modulus, shape)


cdef class RheologyBase(PhysicsBase):
    """Abstract base for rheology models; owns the most-derived C++ model object."""

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "RheologyBase is abstract; instantiate a concrete model "
            "(Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade, Sundberg)."
        )

    def __dealloc__(self):
        self._rheology_ptr.reset()
        self._ptr = NULL

    cdef void _adopt(self, unique_ptr[c_RheologyBase]& model) noexcept:
        """Take ownership of ``model``; the inherited ``_ptr`` observes it."""
        self._rheology_ptr = move(model)
        self._ptr = <c_TidalPyBaseClass*>self._rheology_ptr.get()

    def calc_complex_modulus(self, double modulus,
                             double viscosity,
                             double frequency) -> complex:
        """Complex (shear or bulk) modulus mu* [Pa] at the given forcing frequency.

        Parameters
        ----------
        modulus : float
            Unrelaxed (static) modulus [Pa].
        viscosity : float
            Reference dynamic viscosity [Pa·s].
        frequency : float
            Tidal forcing frequency [rad s-1].

        Returns
        -------
        complex
            Complex modulus [Pa]; real = storage (in-phase), imag = loss (out-of-phase, positive for
            energy loss).

        Notes
        -----
        Assumes a linear viscoelastic regime at a single forcing frequency; the Andrade family further
        assumes a positive forcing frequency.
        """
        self._check_ptr()
        cdef cpp_complex[double] result = self._rheology_ptr.get().calc_complex_modulus(
            modulus, viscosity, frequency)
        return complex(result.real(), result.imag())

    def calc_complex_modulus_vectorize_modulus(self, modulus, viscosity,
                                               double frequency):
        """Complex modulus over equal-length (modulus, viscosity) pairs at one frequency."""
        self._check_ptr()
        return cy_solve_complex_modulus(self._rheology_ptr.get(), modulus, viscosity, frequency, True)

    def calc_complex_modulus_vectorize_frequency(self, double modulus,
                                                 double viscosity, frequency):
        """Complex modulus over a frequency sweep at constant modulus and viscosity."""
        self._check_ptr()
        return cy_solve_complex_modulus(self._rheology_ptr.get(), modulus, viscosity, frequency, True)

    def calc_complex_modulus_vectorize_all(self, modulus, viscosity,
                                           frequency):
        """Complex modulus over element-wise (modulus, viscosity, frequency) triples."""
        self._check_ptr()
        return cy_solve_complex_modulus(self._rheology_ptr.get(), modulus, viscosity, frequency, True)


cdef class Elastic(RheologyBase):
    """Purely elastic: ``mu* = modulus + 0j``. No dissipation, frequency independent."""

    def __init__(self):
        cdef c_RheologyConfig config
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Elastic, config)
        self._adopt(model)


cdef class Viscous(RheologyBase):
    """Purely viscous (Newtonian): ``mu* = i * viscosity * frequency``; the static modulus is unused."""

    def __init__(self):
        cdef c_RheologyConfig config
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Viscous, config)
        self._adopt(model)


cdef class Maxwell(RheologyBase):
    """Standard Maxwell body: ``mu* = 1 / J`` with ``J = 1/modulus - i / (viscosity * frequency)``."""

    def __init__(self):
        cdef c_RheologyConfig config
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Maxwell, config)
        self._adopt(model)


cdef class Voigt(RheologyBase):
    """Voigt-Kelvin element.

    Parameters
    ----------
    voigt_modulus_frac : float, optional
        Voigt element modulus as a multiple of the layer modulus. Default ``5.0``.
    voigt_viscosity_frac : float, optional
        Voigt viscosity as a fraction of the layer viscosity. Default ``0.02``.
    """

    def __init__(self, double voigt_modulus_frac=5.0,
                 double voigt_viscosity_frac=0.02):
        cdef c_RheologyConfig config
        config.voigt_modulus_frac   = voigt_modulus_frac
        config.voigt_viscosity_frac = voigt_viscosity_frac
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Voigt, config)
        self._adopt(model)

    @property
    def voigt_modulus_frac(self) -> float:
        """Voigt modulus fraction [dimensionless] (Voigt modulus as a multiple of the layer modulus)."""
        self._check_ptr()
        return (<c_Voigt*>self._rheology_ptr.get()).get_voigt_modulus_frac()

    @property
    def voigt_viscosity_frac(self) -> float:
        """Voigt viscosity fraction [dimensionless]."""
        self._check_ptr()
        return (<c_Voigt*>self._rheology_ptr.get()).get_voigt_viscosity_frac()


cdef class Burgers(RheologyBase):
    """Burgers rheology: Maxwell and Voigt elements in series.

    Parameters
    ----------
    voigt_modulus_frac : float, optional
        Voigt element modulus as a multiple of the layer modulus. Default ``5.0``.
    voigt_viscosity_frac : float, optional
        Voigt viscosity as a fraction of the layer viscosity. Default ``0.02``.
    """

    def __init__(self, double voigt_modulus_frac=5.0,
                 double voigt_viscosity_frac=0.02):
        cdef c_RheologyConfig config
        config.voigt_modulus_frac   = voigt_modulus_frac
        config.voigt_viscosity_frac = voigt_viscosity_frac
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Burgers, config)
        self._adopt(model)

    @property
    def voigt_modulus_frac(self) -> float:
        """Voigt modulus fraction [dimensionless] (Voigt modulus as a multiple of the layer modulus)."""
        self._check_ptr()
        return (<c_Burgers*>self._rheology_ptr.get()).get_voigt_modulus_frac()

    @property
    def voigt_viscosity_frac(self) -> float:
        """Voigt viscosity fraction [dimensionless]."""
        self._check_ptr()
        return (<c_Burgers*>self._rheology_ptr.get()).get_voigt_viscosity_frac()


cdef class Andrade(RheologyBase):
    """Andrade rheology: a Maxwell body plus a transient term proportional to omega^{-alpha}.

    Parameters
    ----------
    alpha : float, optional
        Andrade exponent [dimensionless]. Default ``0.3``.
    zeta : float, optional
        Andrade timescale ratio [dimensionless]. Default ``1.0``.
    """

    def __init__(self, double alpha=0.3, double zeta=1.0):
        cdef c_RheologyConfig config
        config.alpha = alpha
        config.zeta  = zeta
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Andrade, config)
        self._adopt(model)

    @property
    def alpha(self) -> float:
        """Andrade exponent [dimensionless]."""
        self._check_ptr()
        return (<c_Andrade*>self._rheology_ptr.get()).get_alpha()

    @property
    def zeta(self) -> float:
        """Andrade timescale ratio [dimensionless]."""
        self._check_ptr()
        return (<c_Andrade*>self._rheology_ptr.get()).get_zeta()


cdef class Sundberg(RheologyBase):
    """Sundberg-Cooper rheology: Andrade and Voigt elements summed.

    Parameters
    ----------
    alpha : float, optional
        Andrade exponent [dimensionless]. Default ``0.3``.
    zeta : float, optional
        Andrade timescale ratio [dimensionless]. Default ``1.0``.
    voigt_modulus_frac : float, optional
        Voigt element modulus as a multiple of the layer modulus. Default ``5.0``.
    voigt_viscosity_frac : float, optional
        Voigt viscosity as a fraction of the layer viscosity. Default ``0.02``.
    """

    def __init__(self, double alpha=0.3, double zeta=1.0,
                 double voigt_modulus_frac=5.0,
                 double voigt_viscosity_frac=0.02):
        cdef c_RheologyConfig config
        config.alpha                = alpha
        config.zeta                 = zeta
        config.voigt_modulus_frac   = voigt_modulus_frac
        config.voigt_viscosity_frac = voigt_viscosity_frac
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology(c_RheologyModel.Sundberg, config)
        self._adopt(model)

    @property
    def alpha(self) -> float:
        """Andrade exponent [dimensionless]."""
        self._check_ptr()
        return (<c_Sundberg*>self._rheology_ptr.get()).get_alpha()

    @property
    def zeta(self) -> float:
        """Andrade timescale ratio [dimensionless]."""
        self._check_ptr()
        return (<c_Sundberg*>self._rheology_ptr.get()).get_zeta()

    @property
    def voigt_modulus_frac(self) -> float:
        """Voigt modulus fraction [dimensionless] (Voigt modulus as a multiple of the layer modulus)."""
        self._check_ptr()
        return (<c_Sundberg*>self._rheology_ptr.get()).get_voigt_modulus_frac()

    @property
    def voigt_viscosity_frac(self) -> float:
        """Voigt viscosity fraction [dimensionless]."""
        self._check_ptr()
        return (<c_Sundberg*>self._rheology_ptr.get()).get_voigt_viscosity_frac()


# Every config key any rheology model reads; make_rheology rejects anything else.
RHEOLOGY_CONFIG_KEYS = frozenset({"alpha", "zeta", "voigt_modulus_frac", "voigt_viscosity_frac"})


# The wrapper class of each c_RheologyModel, in enum order.
_RHEOLOGY_CLASSES = (Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade, Sundberg)


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return (
        c_rheology_model_from_name(table_name.encode("utf-8"))
        == c_rheology_model_from_name(model_name.encode("utf-8")))


def make_rheology(str model_name, dict config=None):
    """Build a rheology model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias: ``elastic`` (``off``), ``viscous`` (``newton``), ``voigt``
        (``voigt-kelvin``), ``maxwell``, ``burgers``, ``andrade``, ``sundberg`` (``sundberg-cooper``).
    config : dict, optional
        Model parameters (see ``RHEOLOGY_CONFIG_KEYS``); missing keys fall back to the model defaults.
        ``None`` takes ``[layers.default.shear_rheology]`` from ``TidalPy_Configs_x.toml`` when that table
        names this model, matching what the world builder would attach; an empty dict asks for the
        model's own defaults.

    Returns
    -------
    RheologyBase

    Raises
    ------
    ValueError
        Unknown model name, or a config key that no rheology model reads.
    """
    # None falls back to the same defaults the world-attached path uses.
    config = cy_resolve_factory_config(
        config, "shear_rheology", RHEOLOGY_CONFIG_KEYS, model_name, _same_model, "rheology")

    # The default-constructed config carries the C++ defaults, so only override what the caller gave.
    cdef c_RheologyConfig cfg
    cfg.alpha                = config.get("alpha", cfg.alpha)
    cfg.zeta                 = config.get("zeta", cfg.zeta)
    cfg.voigt_modulus_frac   = config.get("voigt_modulus_frac", cfg.voigt_modulus_frac)
    cfg.voigt_viscosity_frac = config.get("voigt_viscosity_frac", cfg.voigt_viscosity_frac)

    cdef c_RheologyModel model = c_rheology_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_RheologyBase] ptr = c_find_rheology(model, cfg)
    wrapper_class = _RHEOLOGY_CLASSES[<int>model]
    cdef RheologyBase wrapper = wrapper_class.__new__(wrapper_class)
    wrapper._adopt(ptr)
    return wrapper


# Convenience functions. Each builds a stack-allocated C++ model that dies with the call. ``modulus``,
# ``viscosity``, and ``frequency`` accept floats or ndarrays broadcast together; model parameters such
# as ``alpha`` stay scalar.

def elastic(modulus, viscosity, frequency):
    """Complex shear/bulk modulus for the Elastic model [Pa]."""
    cdef c_RheologyConfig cfg
    cdef c_Elastic model = c_Elastic(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)


def viscous(modulus, viscosity, frequency):
    """Complex shear/bulk modulus for the Viscous (Newton) model [Pa]."""
    cdef c_RheologyConfig cfg
    cdef c_Viscous model = c_Viscous(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)


def maxwell(modulus, viscosity, frequency):
    """Complex shear/bulk modulus for the Maxwell model [Pa]."""
    cdef c_RheologyConfig cfg
    cdef c_Maxwell model = c_Maxwell(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)


def voigt(
        modulus,
        viscosity,
        frequency,
        double voigt_modulus_frac=5.0,
        double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Voigt-Kelvin model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.voigt_modulus_frac   = voigt_modulus_frac
    cfg.voigt_viscosity_frac = voigt_viscosity_frac
    cdef c_Voigt model = c_Voigt(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)


def burgers(
        modulus,
        viscosity,
        frequency,
        double voigt_modulus_frac=5.0,
        double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Burgers model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.voigt_modulus_frac   = voigt_modulus_frac
    cfg.voigt_viscosity_frac = voigt_viscosity_frac
    cdef c_Burgers model = c_Burgers(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)


def andrade(
        modulus,
        viscosity,
        frequency,
        double alpha=0.3,
        double zeta=1.0):
    """Complex shear/bulk modulus for the Andrade model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.alpha = alpha
    cfg.zeta  = zeta
    cdef c_Andrade model = c_Andrade(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)


def sundberg(
        modulus,
        viscosity,
        frequency,
        double alpha=0.3,
        double zeta=1.0,
        double voigt_modulus_frac=5.0,
        double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Sundberg-Cooper model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.alpha                = alpha
    cfg.zeta                 = zeta
    cfg.voigt_modulus_frac   = voigt_modulus_frac
    cfg.voigt_viscosity_frac = voigt_viscosity_frac
    cdef c_Sundberg model = c_Sundberg(cfg)
    return cy_solve_complex_modulus(<c_RheologyBase*>&model, modulus, viscosity, frequency, False)
