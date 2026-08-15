# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
classes.pyx
Cython/Python wrappers for TidalPy's base class hierarchy.

TidalPyBaseClass: abstract base wrapper (not directly instantiable from Python)
StructureBase:    spherical geometry base (radius, mass, geometry calcs)
PhysicsBase:      physics model base (model_name, layer observer pointer)
"""

import os as _os

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport make_unique
from libcpp.string cimport string

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address

# Wire this DLL's logger pointer to the shared TidalPy logger so that
# TIDALPY_LOG_* calls inside tidalpy_base_.hpp/binary_.hpp reach the
# correct spdlog instance.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())

# Wire this DLL's config pointer so that tidalpy_config_ptr->d_G and
# other runtime constants resolve to the shared TidalPyConfig instance.
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# TidalPyBaseClass
# =====================================================================================================================
cdef class TidalPyBaseClass:
    """Abstract base for all TidalPy C++ class wrappers.

    Provides binary save/load and schema version access.
    Not directly instantiable — instantiate StructureBase or PhysicsBase instead.
    """

    def __cinit__(self):
        self._ptr = NULL

    cdef void _check_ptr(self) except *:
        if self._ptr is NULL:
            raise RuntimeError(
                f"{type(self).__name__} is not properly initialized."
            )

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
        """Return a dict of this object's configuration.

        Base implementation returns an empty dict.
        Subclasses override to include their own config values.

        Returns
        -------
        dict
        """
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
    """Spherical geometry base class.

    Stores radius [m] and mass [kg]. Provides pure-function geometry
    calculation methods (all const; arguments are explicit, not implicit).

    Parameters
    ----------
    radius_m : float
        Radius in meters (MKS).
    mass_kg : float
        Mass in kilograms (MKS).
    """

    def __cinit__(self, *args, **kwargs):
        # Set _ptr now; address of _struct is stable for the lifetime of this object.
        # Subclasses override _ptr in their own __init__ after constructing a deeper object.
        self._ptr = &self._struct

    def __init__(self, double radius_m, double mass_kg):
        self._struct = c_StructureBase(radius_m, mass_kg)

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
        """Return configuration dict with radius_m and mass_kg [MKS]."""
        return {
            "radius_m": self._struct.get_radius(),
            "mass_kg":  self._struct.get_mass(),
        }


# =====================================================================================================================
# PhysicsBase
# =====================================================================================================================
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
    Subclasses (e.g. rheology models) own the most-derived C++ object
    themselves and leave ``_physics_ptr`` NULL; they set the inherited
    ``_ptr`` to the derived object. The ``model_name`` property therefore
    reads through ``_ptr`` (cast to ``c_PhysicsBase*``) so it works for both
    direct instances and subclasses.
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_PhysicsBase> auto-inits to nullptr; subclasses own their object

    def __init__(self, str model_name):
        cdef string name = model_name.encode("utf-8")
        self._physics_ptr = make_unique[c_PhysicsBase](name)
        self._ptr = <c_TidalPyBaseClass*>self._physics_ptr.get()

    def __dealloc__(self):
        # unique_ptr frees the owned object (if any); just clear the observer ptr.
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
        """Return configuration dict with model_name."""
        self._check_ptr()
        return {
            "model_name": (<c_PhysicsBase*>self._ptr).get_model_name().decode("utf-8"),
        }
