# distutils: language = c++
"""
classes.pxd
Cython declarations for TidalPy's base class hierarchy.

Declares C++ classes and Cython extension types so other modules can cimport
the types and call C-speed methods without going through Python dispatch.

Usage in another extension::

    from TidalPy.Utilities_x.classes_x.classes cimport (
        TidalPyBaseClass, StructureBase, PhysicsBase,
        c_TidalPyBaseClass, c_StructureBase, c_PhysicsBase)

Note: cimporting this pxd brings `bool` (as `cpp_bool`) from libcpp into scope.
Never call bool() as a function in the importing .pyx; use `True if x else False`.
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "tidalpy_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_TidalPyBaseClass:
        string get_schema_version_str() const
        cpp_bool check_schema_compatibility(uint8_t major, uint8_t minor) const
        void save_binary(const string& path) except +
        void load_binary(const string& path, cpp_bool force) except +


cdef extern from "structure_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_StructureBase(c_TidalPyBaseClass):
        c_StructureBase()
        c_StructureBase(double radius_m, double mass_kg) except +
        double get_radius() const
        double get_mass() const
        double calc_surface_area(double radius) const
        double calc_volume_sphere(double radius) const
        double calc_volume_shell(double radius_outer, double radius_inner) const
        double calc_surface_gravity(double mass, double radius) const
        double calc_mean_density(double mass, double volume) const
        double calc_escape_velocity(double mass, double radius) const


cdef extern from "config_entry_.hpp" namespace "tidalpy" nogil:
    # Typed configuration entry reported by a C++ physics model (converted to a dict in Cython).
    cdef enum class c_ConfigEntryKind:
        Double
        Int
        Bool
        String
        DoubleList
        StringList

    cdef cppclass c_ConfigEntry:
        string            key
        c_ConfigEntryKind kind
        double            value_double
        long long         value_int
        cpp_bool          value_bool
        string            value_string
        vector[double]    value_double_list
        vector[string]    value_string_list


cdef extern from "physics_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_PhysicsBase(c_TidalPyBaseClass):
        c_PhysicsBase(const string& model_name) except +
        const string& get_model_name() const
        void set_model_name(const string& name)
        vector[c_ConfigEntry] get_config_entries() const


# Convert a C++ physics model's typed config entries into the dict returned by get_config_dict. Shared by the
# Cython wrappers and by the layer and world writers, which reach attached models through raw pointers.
cdef dict cy_config_entries_to_dict(const vector[c_ConfigEntry]& entries)
cdef dict cy_physics_model_config(const c_PhysicsBase* model_ptr)


# =====================================================================================================================
# Cython wrapper class declarations (for cimport by other extensions)
# =====================================================================================================================
cdef class TidalPyBaseClass:
    cdef c_TidalPyBaseClass* _ptr
    cdef void _check_ptr(self) except *
    cpdef dict get_config_dict(self)


cdef class StructureBase(TidalPyBaseClass):
    cdef c_StructureBase _struct
    cpdef dict get_config_dict(self)


cdef class PhysicsBase(TidalPyBaseClass):
    # Owns the c_PhysicsBase only when PhysicsBase is instantiated directly.
    # Subclasses (e.g. rheology models) leave this empty and own the most-derived
    # object themselves (via their own unique_ptr), setting the inherited
    # TidalPyBaseClass._ptr instead.
    cdef unique_ptr[c_PhysicsBase] _physics_ptr
    cpdef dict get_config_dict(self)
