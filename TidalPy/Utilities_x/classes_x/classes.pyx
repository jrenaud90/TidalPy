# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's base class hierarchy: TidalPyBaseClass, StructureBase, PhysicsBase."""

import difflib
import os as _os

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport make_unique
from libcpp.string cimport string
from libcpp.vector cimport vector

from TidalPy.Utilities_x.logging_x.logger cimport (
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
        """Serialize this object to a TidalPy binary file."""
        self._check_ptr()
        self._ptr.save_binary(path.encode("utf-8"))

    def load_binary(self, str path, cpp_bool force=False):
        """Load this object's state from a TidalPy binary file.

        Parameters
        ----------
        path : str
            Source file path.
        force : bool, optional
            Attempt the load even on a schema version mismatch.
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
        """Save this object's configuration to a TOML file."""
        import toml
        cdef dict config = self.get_config_dict()
        with open(path, 'w', encoding='utf-8') as f:
            toml.dump(config, f)


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


cdef class PhysicsBase(TidalPyBaseClass):
    """Physics model base class.

    Parameters
    ----------
    model_name : str
        Physics model name, e.g. 'maxwell' or 'convection'.

    Notes
    -----
    Subclasses own the most-derived C++ object themselves and leave ``_physics_ptr`` NULL, setting the
    inherited ``_ptr`` instead, so ``model_name`` reads through ``_ptr`` cast to ``c_PhysicsBase*``.
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr auto-inits to nullptr; subclasses own their object

    def __init__(self, str model_name):
        cdef string name = model_name.encode("utf-8")
        self._physics_ptr = make_unique[c_PhysicsBase](name)
        self._ptr = <c_TidalPyBaseClass*>self._physics_ptr.get()

    def __dealloc__(self):
        self._physics_ptr.reset()
        self._ptr = NULL

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


def factory_defaults(str section, accepted_keys, model_name=None, same_model=None) -> dict:
    """The parameters a ``make_*`` factory takes when called with no config: the world builder's defaults.

    The world builder resolves a model's parameters through the layer's own table and then the matching
    model table of ``[layers.default]`` in ``TidalPy_Configs_x.toml`` (``[tides]`` for a tide model). A
    model built directly through its factory shares that second tier, so one file describes both. The
    table's parameters are taken only for the model it names: a parameter of one model must not carry
    over to another that reads the same key differently, as an isotope dataset's reference time would
    into a fixed radiogenic rate.

    Parameters
    ----------
    section : str
        The table below ``[layers.default]``, dotted for nesting (``"material.shear_viscosity"``), or
        ``"tides"`` for the top-level tide table.
    accepted_keys : collection of str
        The keys the family reads (the factory's ``*_CONFIG_KEYS``).
    model_name : str, optional
        The model being built. Given, the table is taken only when its ``model`` names it; a table with no
        ``model`` key, such as ``[tides]``, is taken as is.
    same_model : callable, optional
        ``same_model(table_name, model_name) -> bool``, the family's own alias-aware name test; a
        ``ValueError`` from it counts as a different model. Without it the names are compared as lower-case
        strings, so an alias will not match.

    Returns
    -------
    dict
        The parameters found, empty when the config is not loaded, has no such table, or names another model.
    """
    # Deferred: this module is imported while TidalPy initializes, before config_x exists.
    import TidalPy
    cdef dict config_x = getattr(TidalPy, "config_x", None) or {}
    # `table` stays `object`: walking a dotted section can land on a scalar, which the isinstance check
    # below turns into an empty table. `cdef dict` would raise on the assignment first.
    cdef object table
    cdef object named
    cdef bint matches
    cdef set accepted
    cdef str part
    if section == "tides":
        table = config_x.get("tides", {}) or {}
    else:
        table = (config_x.get("layers", {}) or {}).get("default", {}) or {}
        for part in section.split("."):
            table = table.get(part, {}) or {}
            if not isinstance(table, dict):
                return {}
    named = table.get("model", None)
    if model_name is not None and named is not None:
        if same_model is None:
            matches = str(named).lower() == model_name.lower()
        else:
            try:
                matches = bool(same_model(str(named), model_name))
            except ValueError:
                matches = False
        if not matches:
            return {}
    accepted = set(accepted_keys)
    accepted.discard("model")
    return {key: value for key, value in table.items() if key in accepted}


def check_config_keys(dict config, accepted_keys, str family):
    """Raise ``ValueError`` if a physics-model config holds a key that no model in its family reads.

    Every ``make_*`` factory calls this before building a model, so a misspelled key, most often a missing
    unit suffix, fails loudly instead of silently leaving a parameter at its default.

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
    The check is per family, not per model: a key read by a different model of the family passes, because
    the world builder merges material defaults beneath a user's table.
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
