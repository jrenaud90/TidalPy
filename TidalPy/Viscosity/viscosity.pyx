# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's viscosity models.

Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table, so the wrappers
here only name the models and expose the family's calculation.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_TidalPyBaseClass,
    c_ParamMap,
    c_share_physics,
    cy_param_map,
    cy_resolve_factory_config,
)

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class ViscosityBase(PhysicsBase):
    """Base for viscosity models, which give the pre-melt (solid) viscosity [Pa s] at a temperature and pressure.

    Instantiate a concrete model (``ArrheniusViscosity``, ``ReferenceViscosity``, ``ConstantViscosity``) with its
    parameters as keywords, by argument name or config key, or build one by name with ``make_viscosity``.
    ``get_parameter_info()`` lists every parameter with its unit, default, and bounds.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "ViscosityBase is abstract; instantiate a concrete model "
                "(ArrheniusViscosity, ReferenceViscosity, ConstantViscosity) or call make_viscosity.")
        cdef dict merged = dict(config) if config else {}
        merged.update(parameters)
        cdef c_ParamMap param_map = cy_param_map(merged)
        cdef unique_ptr[c_ViscosityBase] model = c_find_viscosity((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_ViscosityBase](move(model)))

    cdef c_ViscosityBase* _viscosity(self) except NULL:
        self._check_ptr()
        return <c_ViscosityBase*>self._ptr

    def calc_viscosity(self, temperature, pressure=0.0):
        """Dynamic viscosity [Pa s].

        Parameters
        ----------
        temperature : float or np.ndarray
            Temperature [K].
        pressure : float or np.ndarray, optional
            Pressure [Pa]; zero by default.

        Returns
        -------
        float or np.ndarray
            A float when both inputs are floats, otherwise an array of their broadcast shape.
        """
        cdef c_ViscosityBase* model_ptr = self._viscosity()
        cdef vector[vector[double]] inputs
        cdef vector[double] viscosity
        cdef object shape = cy_broadcast_inputs((temperature, pressure), inputs, False)
        if shape is None:
            return model_ptr.calc_viscosity(<double>temperature, <double>pressure)
        with nogil:
            model_ptr.calc_viscosity_vectorize(inputs[0], inputs[1], viscosity)
        return cy_vector_to_ndarray(viscosity, shape)


cdef class ConstantViscosity(ViscosityBase):
    """Viscosity independent of temperature and pressure (alias ``"const"``)."""
    MODEL_NAME = "constant"


cdef class ReferenceViscosity(ViscosityBase):
    """Relative-activation law (alias ``"ref"``): eta = eta_ref exp((E_a / R)(1/T - 1/T_ref) + P V_a / (R T))."""
    MODEL_NAME = "reference"


cdef class ArrheniusViscosity(ViscosityBase):
    """Arrhenius flow law (alias ``"arr"``): eta = A sigma^(1-n) d^m exp((E_a + P V_a)/(R T)), times T when
    ``additional_temp_dependence`` is set."""
    MODEL_NAME = "arrhenius"


# The wrapper class of each model, by canonical name.
_VISCOSITY_CLASSES = {cls.MODEL_NAME: cls for cls in (ArrheniusViscosity, ReferenceViscosity, ConstantViscosity)}


def viscosity_model_names() -> tuple:
    """The canonical names of the viscosity models."""
    cdef vector[string] names = c_viscosity_model_names()
    return tuple(name.decode("utf-8") for name in names)


def canonical_viscosity_name(str model_name) -> str:
    """A viscosity model's canonical name from any of its names or aliases (case-insensitive).

    Raises
    ------
    ValueError
        Unknown model name; the message names the closest one.
    """
    return c_viscosity_canonical_name(model_name.encode("utf-8")).decode("utf-8")


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return canonical_viscosity_name(table_name) == canonical_viscosity_name(model_name)


# Canonical model name -> the config keys that model reads, from its parameter table.
_MODEL_CONFIG_KEYS = {
    name: frozenset(entry["key"] for entry in cls().get_parameter_info()) for name, cls in _VISCOSITY_CLASSES.items()}

# Every config key any viscosity model reads.
VISCOSITY_CONFIG_KEYS = frozenset().union(*_MODEL_CONFIG_KEYS.values())


def viscosity_config_keys(str model_name) -> frozenset:
    """The config keys a viscosity model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _MODEL_CONFIG_KEYS[canonical_viscosity_name(model_name)]


def make_viscosity(str model_name, dict config=None) -> ViscosityBase:
    """Build a viscosity model by name, returning the matching subclass.

    Parameters
    ----------
    model_name : str
        ``"arrhenius"`` (``"arr"``), ``"reference"`` (``"ref"``), or ``"constant"`` (``"const"``).
    config : dict, optional
        Model parameters by config key (see each model's ``get_parameter_info()``). Absent keys take the model's
        defaults. ``None`` takes the shear-viscosity defaults of ``[layers.default]`` in the TidalPy configuration.

    Returns
    -------
    ViscosityBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    """
    cdef str canonical = canonical_viscosity_name(model_name)
    if config is None:
        config = cy_resolve_factory_config(
            config, "material.shear_viscosity", _MODEL_CONFIG_KEYS[canonical], canonical, _same_model, "viscosity")
    return _VISCOSITY_CLASSES[canonical](config)
