# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's viscosity models.

Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table, so the wrappers
here only name the models and expose the family's calculation.
"""

from libcpp.string cimport string
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_ThermoPoint,
    c_share_physics,
    cy_collect_parameters,
    cy_param_map,
    cy_wrap_model,
)
from TidalPy.Utilities.classes.families import ModelFamily

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class ViscosityBase(PhysicsBase):
    """Base for viscosity models, which give the pre-melt (solid) viscosity [Pa s] at a temperature and pressure.

    Instantiate a concrete model (``ArrheniusViscosity``, ``ReferenceViscosity``, ``ConstantViscosity``,
    ``InterpolatedViscosity``) with its parameters positionally (in the order ``get_parameter_info()`` lists them) or
    as keywords (argument names or config keys), or build one by name with ``make_viscosity``.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "ViscosityBase is abstract; instantiate a concrete model "
                "(ArrheniusViscosity, ReferenceViscosity, ConstantViscosity, InterpolatedViscosity) "
                "or call make_viscosity.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_ViscosityBase] model = c_find_viscosity((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_ViscosityBase](move(model)))

    cdef c_ViscosityBase* _viscosity(self) except NULL:
        self._check_ptr()
        return <c_ViscosityBase*>self._ptr

    def calc_viscosity(self, temperature, pressure=0.0, radius=d_NAN):
        """Dynamic viscosity [Pa s].

        Parameters
        ----------
        temperature : float or np.ndarray
            Temperature [K].
        pressure : float or np.ndarray, optional
            Pressure [Pa]; zero by default.
        radius : float or np.ndarray, optional
            Radius [m], read only by ``InterpolatedViscosity``.

        Returns
        -------
        float or np.ndarray
            A float when every input is a float, otherwise an array of their broadcast shape.
        """
        cdef c_ViscosityBase* model_ptr = self._viscosity()
        cdef vector[vector[double]] inputs
        cdef vector[double] viscosity
        cdef c_ThermoPoint point
        cdef object shape = cy_broadcast_inputs((temperature, pressure, radius), inputs, False)
        if shape is None:
            point.pressure    = <double>pressure
            point.temperature = <double>temperature
            point.radius      = <double>radius
            return model_ptr.calc_viscosity(point)
        with nogil:
            model_ptr.calc_viscosity_vectorize(inputs[0], inputs[1], inputs[2], viscosity)
        return cy_vector_to_ndarray(viscosity, shape)


cdef class ConstantViscosity(ViscosityBase):
    """Viscosity independent of temperature and pressure (alias ``"const"``)."""
    MODEL_NAME = "constant"


cdef class ReferenceViscosity(ViscosityBase):
    """Relative-activation law (alias ``"ref"``), anchored at a reference temperature and pressure:
    eta = eta_ref exp((E_a + P V_a)/(R T) - (E_a + P_ref V_a)/(R T_ref))."""
    MODEL_NAME = "reference"


cdef class ArrheniusViscosity(ViscosityBase):
    """Arrhenius flow law (alias ``"arr"``): eta = A sigma^(1-n) d^m exp((E_a + P V_a)/(R T)), times T when
    ``additional_temp_dependence`` is set."""
    MODEL_NAME = "arrhenius"


cdef class InterpolatedViscosity(ViscosityBase):
    """A viscosity profile tabulated in radius (aliases ``"interp"``, ``"interpolated"``): ``radius_m`` and
    ``viscosity_pas`` tables, linear between points and held at the end values beyond them."""
    MODEL_NAME = "interpolate"


cdef class CompositeViscosity(ViscosityBase):
    """Several deformation mechanisms acting in parallel (alias ``"parallel"``): 1 / eta = sum of 1 / eta_i, so the
    weakest dominates. Published ice and olivine flow laws combine diffusion creep, dislocation creep, and
    grain-boundary sliding this way.

    Parameters
    ----------
    *mechanisms : ViscosityBase or dict
        The mechanisms, as viscosity models or as config tables (each with a ``model`` key).
    config : dict, optional
        A config table with a ``mechanisms`` list of tables, as ``get_config_dict()`` returns.
    """
    MODEL_NAME = "composite"

    # Config keys beyond the (empty) parameter table.
    EXTRA_CONFIG_KEYS = frozenset({"mechanisms"})

    def __init__(self, *mechanisms, dict config=None):
        cdef list entries = list(mechanisms)
        cdef object entry
        cdef dict table
        if config:
            unknown = sorted(set(config) - {"model", "mechanisms"})
            if unknown:
                raise ValueError(
                    f"TidalPy: viscosity model 'composite' has no parameter '{unknown[0]}'. Accepted: mechanisms.")
            entries.extend(config.get("mechanisms", []) or [])
        cdef vector[shared_ptr[c_PhysicsBase]] models
        cdef ViscosityBase model
        cdef unique_ptr[c_ViscosityBase] composite
        if not entries:
            composite = c_find_viscosity(b"composite", cy_param_map({}))
        else:
            for entry in entries:
                if isinstance(entry, dict):
                    table = dict(entry)
                    if "model" not in table:
                        raise ValueError("TidalPy: each composite viscosity mechanism table needs a 'model' key.")
                    entry = make_viscosity(str(table.pop("model")), table)
                if not isinstance(entry, ViscosityBase):
                    raise TypeError(
                        f"TidalPy: a composite viscosity mechanism must be a viscosity model or a config table, "
                        f"not {type(entry).__name__}.")
                model = <ViscosityBase>entry
                models.push_back(model._model_sptr)
            composite = c_make_composite_viscosity(models)
        self._set_model(c_share_physics[c_ViscosityBase](move(composite)))

    @property
    def mechanisms(self) -> tuple:
        """The mechanisms, as viscosity models."""
        cdef vector[shared_ptr[c_PhysicsBase]] models = (<c_CompositeViscosity*>self._viscosity()).get_mechanism_models()
        cdef size_t model_i
        return tuple(cy_wrap_model(models[model_i]) for model_i in range(models.size()))


# The family's name lookup: the C++ registry's alias-aware, case-insensitive canonical name.
_FAMILY = ModelFamily(
    "viscosity",
    (ArrheniusViscosity, ReferenceViscosity, ConstantViscosity, InterpolatedViscosity, CompositeViscosity),
    lambda model_name: c_viscosity_canonical_name(model_name.encode("utf-8")).decode("utf-8"))

# Every config key any viscosity model reads.
VISCOSITY_CONFIG_KEYS = _FAMILY.config_keys


def viscosity_model_names() -> tuple:
    """The canonical names of the viscosity models."""
    return _FAMILY.model_names()


def canonical_viscosity_name(str model_name) -> str:
    """A viscosity model's canonical name from any of its names or aliases (case-insensitive).

    Raises
    ------
    ValueError
        Unknown model name; the message names the closest one.
    """
    return _FAMILY.canonical_name(model_name)


def viscosity_config_keys(str model_name) -> frozenset:
    """The config keys a viscosity model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _FAMILY.config_keys_of(model_name)


def make_viscosity(str model_name, dict config=None) -> ViscosityBase:
    """Build a viscosity model by name, returning the matching subclass.

    Parameters
    ----------
    model_name : str
        ``"arrhenius"`` (``"arr"``), ``"reference"`` (``"ref"``), ``"constant"`` (``"const"``), or ``"interpolate"``
        (``"interp"``).
    config : dict, optional
        Model parameters by config key (see each model's ``get_parameter_info()``). Absent keys (all of them for
        ``None``) take the model's defaults.

    Returns
    -------
    ViscosityBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    """
    return _FAMILY.make(model_name, config)
