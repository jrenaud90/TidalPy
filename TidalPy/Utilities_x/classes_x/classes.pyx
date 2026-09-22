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

# Wire this DLL's logger pointer so TIDALPY_LOG_* calls in the C++ headers reach the shared spdlog
# instance.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())

# Wire this DLL's config pointer so tidalpy_config_ptr resolves to the shared TidalPyConfig instance.
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# TidalPyBaseClass
# =====================================================================================================================
cdef class TidalPyBaseClass:
    """Abstract base for all TidalPy C++ class wrappers.

    Provides binary save/load and schema version access.
    Not directly instantiable: instantiate StructureBase or PhysicsBase instead.
    """

    def __cinit__(self):
        self._ptr = NULL

    cdef void _check_ptr(self) except *:
        if self._ptr is NULL:
            raise RuntimeError(
                f"This {type(self).__name__} holds no C++ object: it was never initialized, or it was attached to "
                f"a layer or world, which took ownership of it.")

    def get_schema_version_str(self) -> str:
        """Return the schema version string (e.g. '0.2.0')."""
        self._check_ptr()
        return self._ptr.get_schema_version_str().decode("utf-8")

    def save_binary(self, str path):
        """Serialize this object to a TidalPy binary file.

        Parameters
        ----------
        path : str
            Destination file path.

        Raises
        ------
        IOError
            If the file cannot be opened for writing.
        """
        self._check_ptr()
        self._ptr.save_binary(path.encode("utf-8"))

    def load_binary(self, str path, cpp_bool force=False):
        """Load this object's state from a TidalPy binary file.

        Parameters
        ----------
        path : str
            Source file path.
        force : bool, optional
            If True, attempt to load even on schema version mismatch.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        IOError
            If the file is invalid or schema version is incompatible.
        """
        self._check_ptr()
        if not _os.path.isfile(path):
            raise FileNotFoundError(f"No such file: '{path}'")
        try:
            self._ptr.load_binary(path.encode("utf-8"), force)
        except RuntimeError as exc:
            raise IOError(str(exc)) from exc

    cpdef dict get_config_dict(self):
        """Return a dict of this object's configuration (empty on the base class)."""
        return {}

    def save_config(self, str path):
        """Save this object's configuration to a TOML file.

        Parameters
        ----------
        path : str
            Destination file path (should end in .toml).
        """
        import toml
        config = self.get_config_dict()
        with open(path, 'w', encoding='utf-8') as f:
            toml.dump(config, f)


# =====================================================================================================================
# StructureBase
# =====================================================================================================================
cdef class StructureBase(TidalPyBaseClass):
    """Spherical geometry base class storing radius [m] and mass [kg].

    The geometry methods take explicit arguments instead of reading the stored state.

    Parameters
    ----------
    radius : float
        Radius [m].
    mass : float
        Mass [kg].
    """

    def __cinit__(self, *args, **kwargs):
        # The address of _struct is stable for this object's lifetime; subclasses reset _ptr in __init__.
        self._ptr = &self._struct

    def __init__(self, double radius, double mass):
        self._struct = c_StructureBase(radius, mass)

    def __dealloc__(self):
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def radius(self) -> float:
        """Radius [m]."""
        return self._struct.get_radius()

    @property
    def mass(self) -> float:
        """Mass [kg]."""
        return self._struct.get_mass()

    # ------------------------------------------------------------------------------------------------------------------
    # Geometry calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_surface_area(self, double radius) -> float:
        """Surface area of a sphere [m^2].

        Parameters
        ----------
        radius : float
            Radius [m].
        """
        return self._struct.calc_surface_area(radius)

    def calc_volume_sphere(self, double radius) -> float:
        """Volume of a solid sphere [m^3].

        Parameters
        ----------
        radius : float
            Radius [m].
        """
        return self._struct.calc_volume_sphere(radius)

    def calc_volume_shell(self, double radius_outer, double radius_inner) -> float:
        """Volume of a spherical shell [m^3].

        Parameters
        ----------
        radius_outer : float
            Outer radius [m].
        radius_inner : float
            Inner radius [m].
        """
        return self._struct.calc_volume_shell(radius_outer, radius_inner)

    def calc_surface_gravity(self, double mass, double radius) -> float:
        """Surface gravitational acceleration [m/s^2].

        Parameters
        ----------
        mass : float
            Mass [kg].
        radius : float
            Radius [m].
        """
        return self._struct.calc_surface_gravity(mass, radius)

    def calc_mean_density(self, double mass, double volume) -> float:
        """Mean density [kg/m^3].

        Parameters
        ----------
        mass : float
            Mass [kg].
        volume : float
            Volume [m^3].
        """
        return self._struct.calc_mean_density(mass, volume)

    def calc_escape_velocity(self, double mass, double radius) -> float:
        """Escape velocity [m/s].

        Parameters
        ----------
        mass : float
            Mass [kg].
        radius : float
            Radius [m].
        """
        return self._struct.calc_escape_velocity(mass, radius)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------

    cpdef dict get_config_dict(self):
        """Return configuration dict with radius and mass [MKS]."""
        return {
            "radius_m": self._struct.get_radius(),
            "mass_kg":  self._struct.get_mass(),
        }


# =====================================================================================================================
# PhysicsBase
# =====================================================================================================================
cdef dict cy_config_entries_to_dict(const vector[c_ConfigEntry]& entries):
    """Convert the typed config entries a C++ physics model reports into a Python dict."""
    cdef dict out = {}
    cdef size_t i, j
    cdef str key
    cdef list values
    for i in range(entries.size()):
        key = entries[i].key.decode("utf-8")
        if entries[i].kind == c_ConfigEntryKind.Double:
            out[key] = entries[i].value_double
        elif entries[i].kind == c_ConfigEntryKind.Int:
            out[key] = entries[i].value_int
        elif entries[i].kind == c_ConfigEntryKind.Bool:
            out[key] = True if entries[i].value_bool else False
        elif entries[i].kind == c_ConfigEntryKind.String:
            out[key] = entries[i].value_string.decode("utf-8")
        elif entries[i].kind == c_ConfigEntryKind.DoubleList:
            values = []
            for j in range(entries[i].value_double_list.size()):
                values.append(entries[i].value_double_list[j])
            out[key] = values
        else:
            values = []
            for j in range(entries[i].value_string_list.size()):
                values.append(entries[i].value_string_list[j].decode("utf-8"))
            out[key] = values
    return out


cdef dict cy_physics_model_config(const c_PhysicsBase* model_ptr):
    """Config dict of any C++ physics model (an empty dict for a null pointer).

    Used by the Cython wrappers and by the layer and world writers, which hold their attached models through
    raw pointers.
    """
    cdef vector[c_ConfigEntry] entries
    if model_ptr == NULL:
        return {}
    entries = model_ptr.get_config_entries()
    return cy_config_entries_to_dict(entries)


cdef class PhysicsBase(TidalPyBaseClass):
    """Physics model base class.

    Stores a model name string and a non-owning observer pointer to the layer
    that contains this physics object (set by the layer after construction).

    Parameters
    ----------
    model_name : str
        Physics model name (e.g. 'maxwell', 'convection').

    Notes
    -----
    Subclasses own the most-derived C++ object themselves and leave ``_physics_ptr`` NULL, setting the
    inherited ``_ptr`` instead, so ``model_name`` reads through ``_ptr`` cast to ``c_PhysicsBase*``.
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_PhysicsBase> auto-inits to nullptr; subclasses own their object

    def __init__(self, str model_name):
        cdef string name = model_name.encode("utf-8")
        self._physics_ptr = make_unique[c_PhysicsBase](name)
        self._ptr = <c_TidalPyBaseClass*>self._physics_ptr.get()

    def __dealloc__(self):
        self._physics_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def model_name(self) -> str:
        """Physics model name."""
        self._check_ptr()
        return (<c_PhysicsBase*>self._ptr).get_model_name().decode("utf-8")

    @model_name.setter
    def model_name(self, str value):
        self._check_ptr()
        (<c_PhysicsBase*>self._ptr).set_model_name(value.encode("utf-8"))

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return the configuration dict reported by the C++ model.

        The base entry is the ``model`` name (the key the world builder reads for every physics-model table);
        each C++ subclass appends its own parameters through ``append_config_entries``, so wrapper classes
        never override this method.
        """
        self._check_ptr()
        return cy_physics_model_config(<const c_PhysicsBase*>self._ptr)


# =====================================================================================================================
# Physics-model config key checking
# =====================================================================================================================
def factory_defaults(str section, accepted_keys, model_name=None) -> dict:
    """The parameters a ``make_*`` factory takes when it is called with no config: the world builder's defaults.

    A model attached by the world builder resolves its parameters through the layer's own table and then the
    matching model table of ``[layers.default]`` in ``TidalPy_Configs_x.toml`` (``[tides]`` for a tide model).
    A model built directly through its factory shares that second tier, so the same configuration file
    describes both. The table names one model of the family, and its parameters are taken only for that model
    (compared case-insensitively with ``model_name``): a parameter of one model must not carry over to another
    that reads the same key with a different meaning, as an isotope dataset's reference time would into a
    fixed radiogenic rate. The ``model`` key itself and anything the family does not read are left out.

    Parameters
    ----------
    section : str
        The table below ``[layers.default]``, with dots for nesting (``"shear_rheology"``,
        ``"material.shear_viscosity"``), or ``"tides"`` for the top-level tide table.
    accepted_keys : collection of str
        The keys the family reads (the factory's ``*_CONFIG_KEYS``).
    model_name : str, optional
        The model being built. Given, the table is taken only when its ``model`` names it (a table with no
        ``model`` key, such as ``[tides]``, is taken as is); left out, the table is taken whatever it names.

    Returns
    -------
    dict
        The parameters found, or an empty dict when the configuration is not loaded, has no such table, or the
        table names another model.
    """
    # Deferred: this module is imported while TidalPy initializes, before config_x exists.
    import TidalPy
    config_x = getattr(TidalPy, "config_x", None) or {}
    if section == "tides":
        table = config_x.get("tides", {}) or {}
    else:
        table = (config_x.get("layers", {}) or {}).get("default", {}) or {}
        for part in section.split("."):
            table = table.get(part, {}) or {}
            if not isinstance(table, dict):
                return {}
    named = table.get("model", None)
    if model_name is not None and named is not None and str(named).lower() != model_name.lower():
        return {}
    accepted = set(accepted_keys)
    accepted.discard("model")
    return {key: value for key, value in table.items() if key in accepted}


def check_config_keys(dict config, accepted_keys, str family):
    """Raise ``ValueError`` if a physics-model config holds a key that no model in its family reads.

    Every ``make_*`` factory calls this before building a model, so a misspelled key, most often a missing unit
    suffix, fails loudly instead of silently leaving a parameter at its default.

    Parameters
    ----------
    config : dict or None
        The configuration dict passed to the factory.
    accepted_keys : collection of str
        Every key that some model in the family reads. ``model`` is always accepted as well, so a
        ``get_config_dict()`` result can be passed straight back to its factory.
    family : str
        Model family named in the error message, for example ``"viscosity"``.

    Raises
    ------
    ValueError
        If ``config`` holds any other key. The message names the closest accepted key for each rejected one and
        lists every accepted key.

    Assumptions
    -----------
    The check is per family, not per model: a key read by a different model of the same family passes, because
    the world builder merges material defaults beneath a user's table.
    """
    if not config:
        return
    accepted = set(accepted_keys)
    accepted.add("model")
    rejected = sorted(str(key) for key in config if key not in accepted)
    if not rejected:
        return
    descriptions = []
    for key in rejected:
        close_matches = difflib.get_close_matches(key, accepted, n=1)
        if close_matches:
            descriptions.append(f"'{key}' (did you mean '{close_matches[0]}'?)")
        else:
            descriptions.append(f"'{key}'")
    raise ValueError(
        f"TidalPy: unrecognized {family} config key(s): {', '.join(descriptions)}. "
        f"Accepted keys: {', '.join(sorted(accepted))}.")
