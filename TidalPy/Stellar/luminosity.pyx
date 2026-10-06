# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's stellar luminosity models.

Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table, so the wrappers here
only name the models and expose the family's luminosity relations.

References
----------
- Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 - low-mass mass-luminosity polynomial exponent.
"""

from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp

# The shared vector helpers build their result arrays through the NumPy C API.
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
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
    """Base for stellar luminosity models, which give a star's luminosity from its mass (``calc_luminosity``).

    Instantiate a concrete model (``FixedLuminosity``, ``MassToLuminosity``, ``PowerLawLuminosity``) with its
    parameters positionally (in the order ``get_parameter_info()`` lists them) or as keywords (argument names or config
    keys), or build one by name with ``make_luminosity``. A star shares the model it is given
    (``Star.set_luminosity_model``).
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "LuminosityBase is abstract; instantiate a concrete model (FixedLuminosity, MassToLuminosity, "
                "PowerLawLuminosity) or call make_luminosity.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_LuminosityBase] model = c_find_luminosity((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_LuminosityBase](move(model)))

    cdef c_LuminosityBase* _luminosity(self) except NULL:
        self._check_ptr()
        return <c_LuminosityBase*>self._ptr

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
        return cy_solve_luminosity(self._luminosity(), mass)

    def calc_luminosity_from_temperature(self, double temperature, double radius) -> float:
        """Stefan-Boltzmann luminosity [W] = 4*pi*R^2*sigma*T^4; NaN for non-positive inputs."""
        return self._luminosity().calc_luminosity_from_temperature(temperature, radius)

    def calc_temperature_from_luminosity(self, double luminosity, double radius) -> float:
        """Effective temperature ``T = (L / (4*pi*R^2*sigma))^(1/4)`` [K]; NaN for non-positive inputs."""
        return self._luminosity().calc_temperature_from_luminosity(luminosity, radius)

    def calc_effective_temperature(self, double mass, double radius) -> float:
        """Effective temperature [K] derived from the stellar mass (mass -> L -> T)."""
        return self._luminosity().calc_effective_temperature(mass, radius)


cdef class FixedLuminosity(LuminosityBase):
    """Luminosity supplied directly, the same at every mass (alias ``"constant"``)."""
    MODEL_NAME = "fixed"


cdef class MassToLuminosity(LuminosityBase):
    """Piecewise main-sequence L(M) (aliases ``"cuntz_wang"``, ``"cw"``): the Cuntz and Wang (2018) polynomial
    exponent at low mass, the standard power-law regimes elsewhere. Takes no parameters."""
    MODEL_NAME = "mass_to_luminosity"


cdef class PowerLawLuminosity(LuminosityBase):
    """Single power law ``L = Lsun * coeff * (M / Msun)^exponent`` (alias ``"powerlaw"``)."""
    MODEL_NAME = "power_law"


# The family's name lookup: the C++ registry's alias-aware, case-insensitive canonical name.
_FAMILY = ModelFamily(
    "luminosity",
    (FixedLuminosity, MassToLuminosity, PowerLawLuminosity),
    lambda model_name: c_luminosity_canonical_name(model_name.encode("utf-8")).decode("utf-8"))

# Every config key any luminosity model reads.
LUMINOSITY_CONFIG_KEYS = _FAMILY.config_keys


def luminosity_model_names() -> tuple:
    """The canonical names of the luminosity models."""
    return _FAMILY.model_names()


def luminosity_config_keys(str model_name) -> frozenset:
    """The config keys a luminosity model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _FAMILY.config_keys_of(model_name)


def make_luminosity(str model_name, dict config=None):
    """Build a luminosity model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        ``fixed`` (``constant``), ``mass_to_luminosity`` (``cuntz_wang``, ``cw``), or ``power_law`` (``powerlaw``).
    config : dict, optional
        Model parameters by config key (see each model's ``get_parameter_info()``); absent keys (all of them for
        ``None``) take the model's defaults.

    Returns
    -------
    LuminosityBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    """
    return _FAMILY.make(model_name, config)


# Convenience functions. Each builds the model for the one call; ``mass`` is a float or an ndarray.

cdef object cy_direct_luminosity(str model_name, dict parameters, object mass):
    """Luminosity from a model built for the one call (by name and parameters)."""
    cdef unique_ptr[c_LuminosityBase] model = c_find_luminosity(model_name.encode("utf-8"), cy_param_map(parameters))
    return cy_solve_luminosity(model.get(), mass)


def fixed(mass, double luminosity=0.0):
    """Luminosity for the Fixed model [W] (returns ``luminosity`` regardless of mass)."""
    return cy_direct_luminosity("fixed", {"luminosity": luminosity}, mass)


def mass_to_luminosity(mass):
    """Luminosity for the MassToLuminosity model [W] (piecewise main-sequence relation)."""
    return cy_direct_luminosity("mass_to_luminosity", {}, mass)


def power_law(mass, double coeff=1.0, double exponent=3.5):
    """Luminosity for the PowerLaw model [W]: ``L = Lsun * coeff * (M / Msun)^exponent``."""
    return cy_direct_luminosity("power_law", {"coeff": coeff, "exponent": exponent}, mass)
