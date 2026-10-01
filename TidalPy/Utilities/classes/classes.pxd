# distutils: language = c++
"""Cython declarations for TidalPy's base class hierarchy.

Warning: cimporting this pxd brings ``bool`` (as ``cpp_bool``) from libcpp into scope, so never call
bool() as a function in the importing .pyx; use ``True if x else False``.
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.map cimport map as cpp_map
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t


cdef extern from "tidalpy_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_TidalPyBaseClass:
        string get_schema_version_str() const
        cpp_bool check_schema_compatibility(uint8_t major, uint8_t minor) const
        void save_binary(const string& path) except +
        void load_binary(const string& path, cpp_bool force) except +


cdef extern from "structure_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_StructureBase(c_TidalPyBaseClass):
        c_StructureBase()
        c_StructureBase(double radius, double mass) except +
        double get_radius() const
        double get_mass() const
        double calc_surface_area(double radius) const
        double calc_volume_sphere(double radius) const
        double calc_volume_shell(double radius_outer, double radius_inner) const
        double calc_surface_gravity(double mass, double radius) const
        double calc_mean_density(double mass, double volume) const
        double calc_escape_velocity(double mass, double radius) const


cdef extern from "config_entry_.hpp" namespace "tidalpy" nogil:
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


cdef extern from "param_map_.hpp" namespace "tidalpy" nogil:
    ctypedef cpp_map[string, vector[double]] c_ParamMap

    cdef enum class c_ParamKind:
        Double
        Integer
        Boolean
        Doubles

    cdef enum class c_ParamBounds:
        Any
        Finite
        Positive
        NonNegative
        UnitInterval

    cdef cppclass c_ParamInfo:
        string        name
        string        key
        c_ParamKind   kind
        double        default_value
        vector[double] default_table
        c_ParamBounds bounds
        string        doc


cdef extern from "thermo_point_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_ThermoPoint:
        double pressure
        double temperature
        double radius


cdef extern from "physics_base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_PhysicsBase(c_TidalPyBaseClass):
        c_PhysicsBase(const string& model_name) except +
        const string& get_model_name() const
        vector[c_ConfigEntry] get_config_entries() const
        vector[c_ParamInfo] get_parameter_info() except +
        vector[double] get_parameter(const string& name_or_key) except +
        unique_ptr[c_PhysicsBase] clone_physics() except +
        unique_ptr[c_PhysicsBase] with_parameters(const c_ParamMap& changes) except +

    # A family's unique_ptr as the shared base pointer a wrapper holds; called with move().
    shared_ptr[c_PhysicsBase] c_share_physics[T](unique_ptr[T] model)


# Python parameters (argument names or config keys to floats, booleans, integers, or sequences) as a c_ParamMap.
cdef c_ParamMap cy_param_map(dict parameters) except *

# A spec model's constructor arguments (config, positional in table order, keywords) as one dict.
cdef dict cy_collect_parameters(object model_class, tuple args, dict config, dict parameters)

# Shared by the Cython wrappers and by the layer and world writers, which reach attached models through
# raw pointers.
cdef dict cy_config_entries_to_dict(const vector[c_ConfigEntry]& entries)
cdef dict cy_physics_model_config(const c_PhysicsBase* model_ptr)

# The config a family's make_* factory builds from: the caller's, or the world builder's defaults when it gave none,
# after the family's key check.
cdef dict cy_resolve_factory_config(
    dict config, str section, object accepted_keys, str model_name, object same_model, str family)


cdef class TidalPyBaseClass:
    cdef c_TidalPyBaseClass* _ptr
    cdef void _check_ptr(self) except *
    cpdef dict get_config_dict(self)


cdef class StructureBase(TidalPyBaseClass):
    cdef c_StructureBase _struct
    cpdef dict get_config_dict(self)


cdef class PhysicsBase(TidalPyBaseClass):
    # The model, shared: models are not changed in place (with_parameters returns a new one), so a layer or a solve
    # can hold the same object. A family not yet built on parameter specs leaves this empty, owns its model through
    # its own pointer, and sets the inherited _ptr instead.
    cdef shared_ptr[c_PhysicsBase] _model_sptr
    cdef void _set_model(self, shared_ptr[c_PhysicsBase] model) noexcept
    cpdef dict get_config_dict(self)
