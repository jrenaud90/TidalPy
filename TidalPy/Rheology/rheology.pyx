# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's rheology models. Each returns a complex modulus mu* [Pa].

Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table, so the classes here
only name the models; ``get_parameter_info()`` lists a model's parameters.

References
----------
- Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000
- Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
- Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784
- Nowick and Berry (1972), Anelastic Relaxation in Crystalline Solids (the Zener standard linear solid)
- Kanamori and Anderson (1977), Rev. Geophys., DOI: 10.1029/RG015i001p00105; Wahr and Bergen (1986), GJRAS,
  DOI: 10.1111/j.1365-246X.1986.tb06642.x (the seismic Q model's dispersion)
"""

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector
from cpython.complex cimport PyComplex_FromDoubles

cimport numpy as cnp

# The shared vector helpers build their result arrays through the NumPy C API.
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_complex_vector_to_ndarray
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_share_physics,
    cy_collect_parameters,
    cy_param_map,
)
from TidalPy.Utilities.classes.families import ModelFamily

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
        return PyComplex_FromDoubles(scalar_result.real(), scalar_result.imag())
    with nogil:
        model.calc_complex_modulus_vectorize(inputs[0], inputs[1], inputs[2], complex_modulus)
    return cy_complex_vector_to_ndarray(complex_modulus, shape)


cdef class RheologyBase(PhysicsBase):
    """Base for rheology models, which turn a static modulus [Pa], a viscosity [Pa s], and a forcing frequency
    [rad s-1] into a complex modulus [Pa].

    Instantiate a concrete model with its parameters positionally (in the order ``get_parameter_info()`` lists them)
    or as keywords (argument names or config keys), or build one by name with ``make_rheology``.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "RheologyBase is abstract; instantiate a concrete model "
                "(Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade, Sundberg, Zener, SeismicQ) "
                "or call make_rheology.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_RheologyBase] model = c_find_rheology((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_RheologyBase](move(model)))

    cdef c_RheologyBase* _rheology(self) except NULL:
        self._check_ptr()
        return <c_RheologyBase*>self._ptr

    def calc_complex_modulus(self, modulus, viscosity, frequency):
        """Complex (shear or bulk) modulus mu* [Pa] at the given forcing frequency.

        Parameters
        ----------
        modulus : float or np.ndarray
            Unrelaxed (static) modulus [Pa].
        viscosity : float or np.ndarray
            Reference dynamic viscosity [Pa·s]; for ``SeismicQ``, the quality factor at its reference frequency.
        frequency : float or np.ndarray
            Tidal forcing frequency [rad s-1].

        Returns
        -------
        complex or np.ndarray
            Complex modulus [Pa]; real = storage (in-phase), imag = loss (out-of-phase, positive for energy loss). A
            complex when every input is a float, otherwise an array of their broadcast shape.

        Notes
        -----
        Assumes a linear viscoelastic regime at a single forcing frequency; the Andrade family further assumes a
        positive forcing frequency.
        """
        return cy_solve_complex_modulus(self._rheology(), modulus, viscosity, frequency, False)

    def calc_complex_modulus_vectorize_modulus(self, modulus, viscosity, double frequency):
        """Complex modulus over equal-length (modulus, viscosity) pairs at one frequency."""
        return cy_solve_complex_modulus(self._rheology(), modulus, viscosity, frequency, True)

    def calc_complex_modulus_vectorize_frequency(self, double modulus, double viscosity, frequency):
        """Complex modulus over a frequency sweep at constant modulus and viscosity."""
        return cy_solve_complex_modulus(self._rheology(), modulus, viscosity, frequency, True)

    def calc_complex_modulus_vectorize_all(self, modulus, viscosity, frequency):
        """Complex modulus over element-wise (modulus, viscosity, frequency) triples."""
        return cy_solve_complex_modulus(self._rheology(), modulus, viscosity, frequency, True)


cdef class Elastic(RheologyBase):
    """Purely elastic (alias ``"off"``): ``mu* = modulus + 0j``. No dissipation, frequency independent."""
    MODEL_NAME = "elastic"


cdef class Viscous(RheologyBase):
    """Purely viscous, Newtonian (alias ``"newton"``): ``mu* = i * viscosity * frequency``; the static modulus is
    unused."""
    MODEL_NAME = "viscous"


cdef class Maxwell(RheologyBase):
    """Standard Maxwell body: ``mu* = 1 / J`` with ``J = 1/modulus - i / (viscosity * frequency)``."""
    MODEL_NAME = "maxwell"


cdef class Voigt(RheologyBase):
    """Voigt-Kelvin element (alias ``"voigt-kelvin"``): a spring of ``voigt_modulus_frac * modulus`` in parallel with
    a dashpot of ``voigt_viscosity_frac * viscosity``."""
    MODEL_NAME = "voigt"


cdef class Burgers(RheologyBase):
    """Burgers rheology: Maxwell and Voigt elements in series."""
    MODEL_NAME = "burgers"


cdef class Andrade(RheologyBase):
    """Andrade rheology: a Maxwell body plus a transient term proportional to omega^(-alpha), whose timescale is
    ``zeta`` Maxwell times."""
    MODEL_NAME = "andrade"


cdef class Sundberg(RheologyBase):
    """Sundberg-Cooper rheology (alias ``"sundberg-cooper"``): Andrade and Voigt elements in series."""
    MODEL_NAME = "sundberg"


cdef class Zener(RheologyBase):
    """Zener rheology (standard linear solid; aliases ``"sls"``, ``"standard_linear_solid"``): a relaxed spring in
    parallel with a Maxwell arm.

    mu* = r M + (1 - r) M i omega tau / (1 + i omega tau), with tau = viscosity / ((1 - r) M) and r the
    ``relaxed_modulus_frac``. The response is the modulus M at high frequency and relaxes to r M, not to zero, at low
    frequency. ``r = 0`` is Maxwell and ``r = 1`` is elastic.
    """
    MODEL_NAME = "zener"


cdef class SeismicQ(RheologyBase):
    """Seismic Q (aliases ``"constant_q"``, ``"power_law_q"``): the loss comes from a quality factor measured at a
    reference frequency, not from a viscosity.

    The model reads its viscosity input as that quality factor Q_ref. With s = omega_ref / |omega| and the exponent a,

        Q(omega) = Q_ref s^-a,    Re mu* = M / (1 + D(s) / Q_ref),    Im mu* = Re mu* / Q(omega),

    where M is the modulus at the reference frequency and D(s) = cot(a pi / 2) (s^a - 1), or (2 / pi) ln s at a = 0,
    is the dispersion causality ties to that loss (to first order in 1 / Q). ``a = 0`` keeps one Q at every frequency;
    ``a > 0`` lets it fall toward low frequency as an absorption band does. Zero frequency is unforced: M, no loss.
    """
    MODEL_NAME = "seismic_q"


def _canonical_name(str model_name) -> str:
    return c_rheology_canonical_name(model_name.encode("utf-8")).decode("utf-8")


_FAMILY = ModelFamily(
    "rheology",
    (Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade, Sundberg, Zener, SeismicQ),
    _canonical_name)

# Every config key any rheology model reads.
RHEOLOGY_CONFIG_KEYS = _FAMILY.config_keys


def rheology_model_names() -> tuple:
    """The canonical names of the rheology models."""
    return _FAMILY.model_names()


def canonical_rheology_name(str model_name) -> str:
    """A rheology model's canonical name from any of its names or aliases (case-insensitive).

    Raises
    ------
    ValueError
        Unknown model name; the message names the closest one.
    """
    return _FAMILY.canonical_name(model_name)


def rheology_config_keys(str model_name) -> frozenset:
    """The config keys a rheology model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _FAMILY.config_keys_of(model_name)


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return _FAMILY.same_model(table_name, model_name)


def make_rheology(str model_name, dict config=None):
    """Build a rheology model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias: ``elastic`` (``off``), ``viscous`` (``newton``), ``voigt`` (``voigt-kelvin``),
        ``maxwell``, ``burgers``, ``andrade``, ``sundberg`` (``sundberg-cooper``), ``zener`` (``sls``,
        ``standard_linear_solid``), ``seismic_q`` (``constant_q``, ``power_law_q``).
    config : dict, optional
        Model parameters by config key; missing keys (all of them for ``None``) take the model's defaults.

    Returns
    -------
    RheologyBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    """
    return _FAMILY.make(model_name, config)


# Convenience functions. Each builds the model for the one call. ``modulus``, ``viscosity``, and ``frequency`` accept
# floats or ndarrays broadcast together; model parameters such as ``alpha`` stay scalar.

cdef object cy_direct_complex_modulus(str model_name, dict parameters, object modulus, object viscosity,
                                      object frequency):
    cdef unique_ptr[c_RheologyBase] model = c_find_rheology(model_name.encode("utf-8"), cy_param_map(parameters))
    return cy_solve_complex_modulus(model.get(), modulus, viscosity, frequency, False)


def elastic(modulus, viscosity, frequency):
    """Complex shear/bulk modulus for the Elastic model [Pa]."""
    return cy_direct_complex_modulus("elastic", {}, modulus, viscosity, frequency)


def viscous(modulus, viscosity, frequency):
    """Complex shear/bulk modulus for the Viscous (Newton) model [Pa]."""
    return cy_direct_complex_modulus("viscous", {}, modulus, viscosity, frequency)


def maxwell(modulus, viscosity, frequency):
    """Complex shear/bulk modulus for the Maxwell model [Pa]."""
    return cy_direct_complex_modulus("maxwell", {}, modulus, viscosity, frequency)


def voigt(
        modulus,
        viscosity,
        frequency,
        double voigt_modulus_frac=5.0,
        double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Voigt-Kelvin model [Pa]."""
    return cy_direct_complex_modulus(
        "voigt",
        {"voigt_modulus_frac": voigt_modulus_frac, "voigt_viscosity_frac": voigt_viscosity_frac},
        modulus, viscosity, frequency)


def burgers(
        modulus,
        viscosity,
        frequency,
        double voigt_modulus_frac=5.0,
        double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Burgers model [Pa]."""
    return cy_direct_complex_modulus(
        "burgers",
        {"voigt_modulus_frac": voigt_modulus_frac, "voigt_viscosity_frac": voigt_viscosity_frac},
        modulus, viscosity, frequency)


def andrade(
        modulus,
        viscosity,
        frequency,
        double alpha=0.3,
        double zeta=1.0):
    """Complex shear/bulk modulus for the Andrade model [Pa]."""
    return cy_direct_complex_modulus("andrade", {"alpha": alpha, "zeta": zeta}, modulus, viscosity, frequency)


def sundberg(
        modulus,
        viscosity,
        frequency,
        double alpha=0.3,
        double zeta=1.0,
        double voigt_modulus_frac=5.0,
        double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Sundberg-Cooper model [Pa]."""
    return cy_direct_complex_modulus(
        "sundberg",
        {"alpha": alpha, "zeta": zeta, "voigt_modulus_frac": voigt_modulus_frac,
         "voigt_viscosity_frac": voigt_viscosity_frac},
        modulus, viscosity, frequency)


def zener(
        modulus,
        viscosity,
        frequency,
        double relaxed_modulus_frac=0.5):
    """Complex shear/bulk modulus for the Zener (standard linear solid) model [Pa]."""
    return cy_direct_complex_modulus(
        "zener", {"relaxed_modulus_frac": relaxed_modulus_frac}, modulus, viscosity, frequency)


def seismic_q(
        modulus,
        quality_factor,
        frequency,
        double reference_frequency_rad_s=6.283185307179586,
        double q_frequency_exponent=0.0):
    """Complex shear/bulk modulus for the seismic Q model [Pa]; ``quality_factor`` is Q at the reference frequency."""
    return cy_direct_complex_modulus(
        "seismic_q",
        {"reference_frequency_rad_s": reference_frequency_rad_s, "q_frequency_exponent": q_frequency_exponent},
        modulus, quality_factor, frequency)
