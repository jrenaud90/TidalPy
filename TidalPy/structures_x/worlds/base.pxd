# distutils: language = c++
"""
base.pxd
Cython declarations for TidalPy's base world class (Phase 8b).

Exports c_WorldConfig, c_BaseWorld, and the Python wrapper BaseWorld so other
extensions can cimport and use C-speed access.
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr

from TidalPy.Utilities_x.classes_x.classes cimport (
    TidalPyBaseClass,
    StructureBase,
    c_StructureBase,
)


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_WorldConfig:
        string   name
        string   world_type_str
        double   radius_m
        double   mass_kg
        double   albedo
        double   emissivity
        double   obliquity_rad
        double   spin_frequency_rad_s

    cdef cppclass c_BaseWorld(c_StructureBase):
        c_BaseWorld()
        c_BaseWorld(const c_WorldConfig& cfg) except +
        const string& get_name()         const
        const string& get_world_type()   const
        double   get_albedo()            const
        double   get_emissivity()        const
        double   get_obliquity()         const
        double   get_spin_frequency()    const
        double   calc_surface_gravity()  const
        double   calc_escape_velocity()  const
        double   calc_mean_density()     const
        double   calc_equilibrium_temperature(double insolation_flux_w_m2) const
        void     set_spin_frequency(double freq_rad_s)
        void     set_obliquity(double obliq_rad)


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class BaseWorld(StructureBase):
    cdef unique_ptr[c_BaseWorld] _world_ptr   # owns the most-derived C++ world object
    cdef public dict source_config            # normalized config the world was built from (or None)
    cpdef dict get_config_dict(self)
