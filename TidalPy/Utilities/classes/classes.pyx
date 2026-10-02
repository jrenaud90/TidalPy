# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's base class hierarchy (TidalPyBaseClass, StructureBase, PhysicsBase) and the
parameter conversions every spec-driven physics model shares."""

import difflib
import os as _os

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport make_unique, shared_ptr, unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address

# Wire this DLL's logger pointer so TIDALPY_LOG_* calls in the C++ headers reach the shared spdlog.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())

# Wire this DLL's config pointer so tidalpy_config_ptr resolves to the shared TidalPyConfig instance.
set_tidalpy_config_ptr(get_shared_config_address())


cdef class TidalPyBaseClass:
    """Abstract base for all TidalPy C++ class wrappers: binary save/load and schema version access."""

    def __cinit__(self):
        self._ptr = NULL

    cdef void _check_ptr(self) except *:
        if self._ptr is NULL:
            raise RuntimeError(
                f"This {type(self).__name__} holds no C++ object: it was never initialized, or it was attached to "
                f"a layer or world, which took ownership of it.")

    def get_schema_version_str(self) -> str:
        """Schema version string, e.g. '0.2.0'."""
        self._check_ptr()
        return self._ptr.get_schema_version_str().decode("utf-8")

    def save_binary(self, str path):
        """Serialize this object to a TidalPy binary file.

        The record is written to a temporary file beside ``path`` and renamed over ``path`` only once it is complete,
        so a failed save leaves any previous file at ``path`` unchanged.

        Parameters
        ----------
        path : str
            Destination file path.

        Raises
        ------
        IOError
            The file cannot be written, or the finished file cannot replace an existing one at ``path`` (for example
            while another program holds it open on Windows).
        """
        self._check_ptr()
        try:
            self._ptr.save_binary(path.encode("utf-8"))
        except RuntimeError as exc:
            raise IOError(str(exc)) from exc

    def load_binary(self, str path, cpp_bool force=False):
        """Load this object's state from a TidalPy binary file.

        A load that raises leaves every setting the object saves as it was before the call.

        Parameters
        ----------
        path : str
            Source file path.
        force : bool, optional
            Attempt the load even on a schema version mismatch.

        Raises
        ------
        FileNotFoundError
            ``path`` does not exist.
        IOError
            The file holds a record of another class, has an incompatible schema version, was written in another byte
            order, or is corrupt: it ends inside a record, a record size disagrees with what this build reads, or
            bytes are left over after the record.
        """
        self._check_ptr()
        if not _os.path.isfile(path):
            raise FileNotFoundError(f"No such file: '{path}'")
        try:
            self._ptr.load_binary(path.encode("utf-8"), force)
        except RuntimeError as exc:
            raise IOError(str(exc)) from exc

    cpdef dict get_config_dict(self):
        """This object's configuration; empty on the base class."""
        return {}

    def save_config(self, str path):
        """Save this object's configuration (:meth:`get_config_dict`) to a TOML file.

        The file starts with the comment header naming the TidalPy, SciPy, and CyRK versions, as every saved
        configuration does, and is written with LF newlines.
        """
        # Deferred like the TidalPy import in factory_defaults: this module is imported while TidalPy initializes.
        from TidalPy.configurations import write_config_toml
        write_config_toml(self.get_config_dict(), path, f"TidalPy {type(self).__name__} configuration")


cdef class StructureBase(TidalPyBaseClass):
    """Spherical geometry base class storing radius [m] and mass [kg].

    The geometry methods take explicit arguments rather than reading the stored state.

    Parameters
    ----------
    radius : float
        Radius [m].
    mass : float
        Mass [kg].
    """

    def __cinit__(self, *args, **kwargs):
        # _struct's address is stable for this object's lifetime; subclasses reset _ptr in __init__.
        self._ptr = &self._struct

    def __init__(self, double radius, double mass):
        self._struct = c_StructureBase(radius, mass)

    def __dealloc__(self):
        self._ptr = NULL

    @property
    def radius(self) -> float:
        """Radius [m]."""
        return self._struct.get_radius()

    @property
    def mass(self) -> float:
        """Mass [kg]."""
        return self._struct.get_mass()

    def calc_surface_area(self, double radius) -> float:
        """Surface area of a sphere [m^2]."""
        return self._struct.calc_surface_area(radius)

    def calc_volume_sphere(self, double radius) -> float:
        """Volume of a solid sphere [m^3]."""
        return self._struct.calc_volume_sphere(radius)

    def calc_volume_shell(self, double radius_outer, double radius_inner) -> float:
        """Volume of a spherical shell [m^3]."""
        return self._struct.calc_volume_shell(radius_outer, radius_inner)

    def calc_surface_gravity(self, double mass, double radius) -> float:
        """Surface gravitational acceleration [m/s^2]."""
        return self._struct.calc_surface_gravity(mass, radius)

    def calc_mean_density(self, double mass, double volume) -> float:
        """Mean density [kg/m^3]."""
        return self._struct.calc_mean_density(mass, volume)

    def calc_escape_velocity(self, double mass, double radius) -> float:
        """Escape velocity [m/s]."""
        return self._struct.calc_escape_velocity(mass, radius)

    cpdef dict get_config_dict(self):
        """Radius and mass [MKS]."""
        return {
            "radius_m": self._struct.get_radius(),
            "mass_kg":  self._struct.get_mass(),
        }


cdef dict cy_config_entries_to_dict(const vector[c_ConfigEntry]& entries):
    cdef dict out = {}
    cdef size_t i, j
    cdef str key
    cdef list values
    cdef const c_ConfigEntry* entry_ptr = NULL
    for i in range(entries.size()):
        entry_ptr = &entries[i]
        key = entry_ptr.key.decode("utf-8")
        if entry_ptr.kind == c_ConfigEntryKind.Double:
            out[key] = entry_ptr.value_double
        elif entry_ptr.kind == c_ConfigEntryKind.Int:
            out[key] = entry_ptr.value_int
        elif entry_ptr.kind == c_ConfigEntryKind.Bool:
            out[key] = True if entry_ptr.value_bool else False
        elif entry_ptr.kind == c_ConfigEntryKind.String:
            out[key] = entry_ptr.value_string.decode("utf-8")
        elif entry_ptr.kind == c_ConfigEntryKind.DoubleList:
            values = []
            for j in range(entry_ptr.value_double_list.size()):
                values.append(entry_ptr.value_double_list[j])
            out[key] = values
        elif entry_ptr.kind == c_ConfigEntryKind.Table:
            out[key] = cy_config_entries_to_dict(entry_ptr.value_table)
        elif entry_ptr.kind == c_ConfigEntryKind.TableList:
            values = []
            for j in range(entry_ptr.value_table_list.size()):
                values.append(cy_config_entries_to_dict(entry_ptr.value_table_list[j]))
            out[key] = values
        else:
            values = []
            for j in range(entries[i].value_string_list.size()):
                values.append(entries[i].value_string_list[j].decode("utf-8"))
            out[key] = values
    return out


cdef dict cy_physics_model_config(const c_PhysicsBase* model_ptr):
    """Config dict of any C++ physics model; empty for a null pointer."""
    cdef vector[c_ConfigEntry] entries
    if model_ptr == NULL:
        return {}
    entries = model_ptr.get_config_entries()
    return cy_config_entries_to_dict(entries)


cdef c_ParamMap cy_param_map(dict parameters) except *:
    """Python parameters as a c_ParamMap: each value a float, an int, a bool, or a sequence of floats (a table).

    Keys may be argument names or config keys; the C++ model accepts either. ``None`` values are skipped, so a caller
    can pass optional arguments straight through, and so is the ``model`` key, which names the model rather than a
    parameter, so a ``get_config_dict()`` result goes straight back to its factory.
    """
    cdef c_ParamMap param_map
    cdef vector[double] values
    cdef object key
    cdef object value
    cdef object item
    for key, value in parameters.items():
        if value is None or key == "model":
            continue
        if not isinstance(key, str):
            raise TypeError(f"TidalPy: parameter names must be strings; got {key!r}.")
        if isinstance(value, str):
            raise TypeError(f"TidalPy: parameter '{key}' takes a number, not the string {value!r}.")
        values.clear()
        if isinstance(value, (bool, int, float)):
            values.push_back(<double>float(value))
        else:
            try:
                for item in value:
                    values.push_back(<double>float(item))
            except TypeError:
                # A NumPy scalar, or another object that converts to float.
                values.clear()
                values.push_back(<double>float(value))
        param_map[(<str>key).encode("utf-8")] = values
    return param_map


cdef object cy_wrap_model(shared_ptr[c_PhysicsBase] model):
    """Any spec model as an instance of its family's Python class for its model name.

    The class comes from the registry ModelFamily fills (TidalPy.Utilities.classes.families), keyed by the family
    name the C++ model reports. Raises TypeError for a family or model no Python class is registered for.
    """
    if model.get() == NULL:
        return None
    # Deferred: families imports this module.
    from TidalPy.Utilities.classes.families import model_class
    cls = model_class(model.get().get_family_name().decode("utf-8"), model.get().get_model_name().decode("utf-8"))
    cdef PhysicsBase wrapper = cls.__new__(cls)
    wrapper._set_model(model)
    return wrapper


# Spec-model class -> {argument name: config key} of its parameters, read once from a default instance.
_PARAMETER_SPELLINGS = {}


def canonical_parameter_keys(object model_class, dict table) -> dict:
    """``table`` with each parameter of ``model_class`` given by its argument name renamed to its config key.

    Every spec model takes a parameter by either spelling. Merging two tables (a config and keywords, a preset and its
    overrides) by key works only when both use one spelling, so callers canonicalize first. Keys that are not
    parameters of the class (slot tables, ``model``) pass through unchanged.

    Raises
    ------
    TypeError
        The table gives one parameter under both spellings.
    """
    # An empty table needs no spellings; building the default instance that lists them would recurse for a class whose
    # constructor canonicalizes its own (empty) keywords.
    if not table:
        return {}
    if model_class not in _PARAMETER_SPELLINGS:
        _PARAMETER_SPELLINGS[model_class] = {
            entry["name"]: entry["key"] for entry in model_class().get_parameter_info()}
    cdef dict spellings = _PARAMETER_SPELLINGS[model_class]
    cdef dict out = {}
    for key, value in table.items():
        canonical = spellings.get(key, key)
        if canonical in out:
            raise TypeError(f"TidalPy: '{key}' and '{canonical}' are the same parameter of {model_class.__name__}; "
                            "give it once.")
        out[canonical] = value
    return out


cdef dict cy_collect_parameters(object model_class, tuple args, dict config, dict parameters):
    """A spec model's constructor arguments as one dict keyed by config key: ``config``, then the positional ``args`` in
    the order of the model's parameter table, then the keywords (argument names or config keys), each overriding what
    came before.

    Raises
    ------
    TypeError
        More positional arguments than parameters, a parameter given both positionally and by keyword, or one
        parameter given under both of its spellings in one table.
    """
    cdef dict merged = canonical_parameter_keys(model_class, config) if config else {}
    cdef dict keywords = canonical_parameter_keys(model_class, parameters) if parameters else {}
    cdef list entries
    cdef Py_ssize_t arg_i
    if args:
        # A default instance lists the parameters in table order.
        entries = model_class().get_parameter_info()
        if len(args) > len(entries):
            raise TypeError(
                f"{model_class.__name__} takes at most {len(entries)} positional parameters "
                f"({', '.join(entry['name'] for entry in entries) if entries else 'none'}); got {len(args)}.")
        for arg_i in range(len(args)):
            if entries[arg_i]["key"] in keywords:
                raise TypeError(f"{model_class.__name__} got multiple values for '{entries[arg_i]['name']}'.")
            merged[entries[arg_i]["key"]] = args[arg_i]
    merged.update(keywords)
    return merged


cdef object cy_param_value(c_ParamKind kind, const vector[double]& values):
    """A parameter's values as the Python type its kind names."""
    cdef size_t value_i
    if kind == c_ParamKind.Doubles:
        return [values[value_i] for value_i in range(values.size())]
    if values.size() == 0:
        return None
    if kind == c_ParamKind.Boolean:
        return True if values[0] != 0.0 else False
    if kind == c_ParamKind.Integer:
        return int(values[0])
    return values[0]


cdef str cy_param_kind_name(c_ParamKind kind):
    if kind == c_ParamKind.Integer:
        return "int"
    if kind == c_ParamKind.Boolean:
        return "bool"
    if kind == c_ParamKind.Doubles:
        return "list[float]"
    return "float"


cdef str cy_param_bounds_name(c_ParamBounds bounds):
    if bounds == c_ParamBounds.Finite:
        return "finite"
    if bounds == c_ParamBounds.Positive:
        return "positive"
    if bounds == c_ParamBounds.NonNegative:
        return "non-negative"
    if bounds == c_ParamBounds.UnitInterval:
        return "unit interval"
    if bounds == c_ParamBounds.PositiveOrInfinite:
        return "positive or infinite"
    return "any"


cdef class PhysicsBase(TidalPyBaseClass):
    """Physics model base class.

    A model built on a parameter spec answers a generic parameter interface here: ``parameters``,
    ``get_parameter``, ``get_parameter_info``, ``with_parameters``, and attribute access by parameter name
    (``model.reference_viscosity``). Models are not changed in place: ``with_parameters`` returns a new model, so one
    model can be shared by several layers.

    Parameters
    ----------
    model_name : str
        Physics model name, e.g. 'maxwell' or 'convection'.

    Notes
    -----
    A family not yet built on parameter specs owns its model through its own pointer and sets the inherited
    ``_ptr``; the generic parameter methods then report no parameters or raise.
    """

    def __init__(self, str model_name):
        cdef string name = model_name.encode("utf-8")
        cdef unique_ptr[c_PhysicsBase] model = make_unique[c_PhysicsBase](name)
        self._set_model(c_share_physics[c_PhysicsBase](move(model)))

    def __dealloc__(self):
        self._model_sptr.reset()
        self._ptr = NULL

    cdef void _set_model(self, shared_ptr[c_PhysicsBase] model) noexcept:
        """Hold ``model``; the inherited ``_ptr`` observes it."""
        self._model_sptr = model
        self._ptr = <c_TidalPyBaseClass*>self._model_sptr.get()

    def get_parameter_info(self) -> list:
        """Descriptions of the model's parameters, in order.

        Returns
        -------
        list of dict
            One dict per parameter: ``name`` (the argument name), ``key`` (the config key, unit suffix included),
            ``kind`` (``"float"``, ``"int"``, ``"bool"``, or ``"list[float]"``), ``default``, ``bounds``, and
            ``doc``. Empty for a model without a parameter spec.
        """
        self._check_ptr()
        cdef vector[c_ParamInfo] info = (<c_PhysicsBase*>self._ptr).get_parameter_info()
        cdef vector[double] default_values
        cdef list out = []
        cdef size_t info_i
        for info_i in range(info.size()):
            default_values.assign(1, info[info_i].default_value)
            out.append({
                "name":    info[info_i].name.decode("utf-8"),
                "key":     info[info_i].key.decode("utf-8"),
                "kind":    cy_param_kind_name(info[info_i].kind),
                "default": cy_param_value(
                    info[info_i].kind,
                    info[info_i].default_table if info[info_i].kind == c_ParamKind.Doubles else default_values),
                "bounds":  cy_param_bounds_name(info[info_i].bounds),
                "doc":     info[info_i].doc.decode("utf-8"),
            })
        return out

    def get_parameter(self, str name):
        """A parameter's value by argument name or config key; a list for a table.

        Raises
        ------
        ValueError
            The model has no such parameter; the message names the closest one.
        """
        self._check_ptr()
        cdef string encoded = name.encode("utf-8")
        cdef vector[double] values = (<c_PhysicsBase*>self._ptr).get_parameter(encoded)
        cdef vector[c_ParamInfo] info = (<c_PhysicsBase*>self._ptr).get_parameter_info()
        cdef size_t info_i
        for info_i in range(info.size()):
            if info[info_i].name == encoded or info[info_i].key == encoded:
                return cy_param_value(info[info_i].kind, values)
        return cy_param_value(c_ParamKind.Double, values)

    @property
    def parameters(self) -> dict:
        """Every parameter by argument name (no unit suffix); ``get_config_dict`` gives them by config key."""
        self._check_ptr()
        cdef vector[c_ParamInfo] info = (<c_PhysicsBase*>self._ptr).get_parameter_info()
        cdef dict out = {}
        cdef size_t info_i
        for info_i in range(info.size()):
            out[info[info_i].name.decode("utf-8")] = cy_param_value(
                info[info_i].kind, (<c_PhysicsBase*>self._ptr).get_parameter(info[info_i].key))
        return out

    def with_parameters(self, **changes):
        """A new model of the same kind with the given parameters changed (argument names or config keys).

        This model is unchanged, and the new one is validated like a newly built one.

        Raises
        ------
        ValueError
            An unknown parameter (the message names the closest one) or a value outside its bounds.
        """
        self._check_ptr()
        cdef c_ParamMap param_map = cy_param_map(changes)
        cdef unique_ptr[c_PhysicsBase] changed = (<c_PhysicsBase*>self._ptr).with_parameters(param_map)
        cdef PhysicsBase wrapper = type(self).__new__(type(self))
        wrapper._set_model(c_share_physics[c_PhysicsBase](move(changed)))
        return wrapper

    def __copy__(self):
        # Models are not changed in place, so a copy can share this one.
        return self

    def __deepcopy__(self, memo):
        self._check_ptr()
        cdef unique_ptr[c_PhysicsBase] copied = (<c_PhysicsBase*>self._ptr).clone_physics()
        cdef PhysicsBase wrapper = type(self).__new__(type(self))
        wrapper._set_model(c_share_physics[c_PhysicsBase](move(copied)))
        return wrapper

    def load_binary(self, str path, cpp_bool force=False):
        """Load this object's state from a TidalPy binary file.

        A spec model is read into a fresh copy that then replaces this wrapper's model, so a layer that shares the
        previous model keeps it. A load that raises leaves this object as it was.

        Raises
        ------
        FileNotFoundError
            ``path`` does not exist.
        IOError
            The file holds a record of another class, has an incompatible schema version, or is corrupt.
        """
        self._check_ptr()
        cdef unique_ptr[c_PhysicsBase] fresh
        if self._model_sptr.get() != NULL:
            try:
                fresh = (<c_PhysicsBase*>self._ptr).clone_physics()
            except RuntimeError:
                # A model without a parameter spec cannot be copied, so it is read in place.
                pass
        if fresh.get() == NULL:
            return TidalPyBaseClass.load_binary(self, path, force)
        if not _os.path.isfile(path):
            raise FileNotFoundError(f"No such file: '{path}'")
        try:
            fresh.get().load_binary(path.encode("utf-8"), force)
        except RuntimeError as exc:
            raise IOError(str(exc)) from exc
        self._set_model(c_share_physics[c_PhysicsBase](move(fresh)))

    def __getattr__(self, str name):
        # Parameters read as attributes, by argument name or config key. Python calls this only after normal lookup
        # fails, so methods and properties win.
        if name.startswith("_") or self._ptr is NULL:
            raise AttributeError(name)
        cdef vector[c_ParamInfo] info = (<c_PhysicsBase*>self._ptr).get_parameter_info()
        cdef string encoded = name.encode("utf-8")
        cdef size_t info_i
        for info_i in range(info.size()):
            if info[info_i].name == encoded or info[info_i].key == encoded:
                return cy_param_value(
                    info[info_i].kind, (<c_PhysicsBase*>self._ptr).get_parameter(info[info_i].key))
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __dir__(self):
        cdef list names = list(object.__dir__(self))
        if self._ptr is not NULL:
            names.extend(entry["name"] for entry in self.get_parameter_info())
        return names

    @property
    def model_name(self) -> str:
        """Physics model name. Read-only: the name is what the model computes, so a different model is a new
        object built by the family's ``make_*`` factory."""
        self._check_ptr()
        return (<c_PhysicsBase*>self._ptr).get_model_name().decode("utf-8")

    cpdef dict get_config_dict(self):
        """The configuration dict the C++ model reports.

        The base entry is the ``model`` name, the key the world builder reads for every physics-model
        table; each C++ subclass appends its own parameters, so wrapper classes never override this.
        """
        self._check_ptr()
        return cy_physics_model_config(<const c_PhysicsBase*>self._ptr)


def factory_defaults(str section, accepted_keys) -> dict:
    """The keys of the ``[section]`` table of ``TidalPy_Configs.toml`` that a family's ``make_*`` factory reads.

    Two families have such a table: tide models take ``[tides]`` (the per-degree lists the world builder also falls
    back on) under the keys a call gives, and an isotope radiogenics model given no dataset takes the ``isotopes``
    of ``[radiogenics]``. Every other family takes its models' own defaults, which their parameter tables state.

    Parameters
    ----------
    section : str
        The top-level table, ``"tides"`` or ``"radiogenics"``.
    accepted_keys : collection of str
        The keys the family reads (the factory's ``*_CONFIG_KEYS``); the table's other keys are left out.

    Returns
    -------
    dict
        The parameters found, empty when the config is not loaded or has no such table.
    """
    # Deferred: this module is imported while TidalPy initializes, before config exists.
    import TidalPy
    cdef dict config = getattr(TidalPy, "config", None) or {}
    cdef object table = config.get(section, {}) or {}
    if not isinstance(table, dict):
        return {}
    return {key: value for key, value in table.items() if key in accepted_keys and key != "model"}


def check_config_keys(dict config, accepted_keys, str family):
    """Raise ``ValueError`` if a physics-model config holds a key that no model in its family reads.

    The radiogenics, tide, and luminosity factories call this before building a model (the spec-model families
    check their keys in C++), so a misspelled key, most often a missing unit suffix, fails loudly instead of
    silently leaving a parameter at its default.

    Parameters
    ----------
    config : dict or None
        The configuration dict passed to the factory.
    accepted_keys : collection of str
        Every key some model in the family reads. ``model`` is always accepted too, so a
        ``get_config_dict()`` result can go straight back to its factory.
    family : str
        Model family named in the error message, for example ``"viscosity"``.

    Raises
    ------
    ValueError
        The message names the closest accepted key for each rejected one.

    Notes
    -----
    The check is per family, not per model: a key read by a different model of the family passes, so one
    table can carry the settings of several models of the family.
    """
    cdef set accepted
    cdef list rejected
    cdef list close_matches
    cdef str key
    if not config:
        return
    accepted = set(accepted_keys)
    accepted.add("model")
    rejected = sorted(str(given) for given in config if given not in accepted)
    if not rejected:
        return
    cdef list descriptions = []
    for key in rejected:
        close_matches = difflib.get_close_matches(key, accepted, n=1)
        if close_matches:
            descriptions.append(f"'{key}' (did you mean '{close_matches[0]}'?)")
        else:
            descriptions.append(f"'{key}'")
    raise ValueError(
        f"TidalPy: unrecognized {family} config key(s): {', '.join(descriptions)}. "
        f"Accepted keys: {', '.join(sorted(accepted))}.")
