# distutils: language = c++
"""
stellar.pxd
Cython declarations for TidalPy's star world class (Phase 8b).
"""

from libcpp.string cimport string

from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "stellar_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_StarConfig:
        # Inherited from c_WorldConfig:
        string   name
        string   world_type_str
        double   radius_m
        double   mass_kg
        double   albedo
        double   emissivity
        double   obliquity_rad
        double   spin_frequency_rad_s
        # StarWorld additions:
        double   effective_temperature_k
        double   luminosity_w

    cdef cppclass c_StarWorld(c_BaseWorld):
        c_StarWorld()
        c_StarWorld(const c_StarConfig& cfg) except +
        double get_effective_temperature() const
        double get_luminosity()            const
        double calc_luminosity_from_temperature(double temperature_k) const
        double calc_temperature_from_luminosity(double luminosity_w)  const
        void   set_effective_temperature(double temperature_k)
        void   set_luminosity(double luminosity_w)


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class StarWorld(BaseWorld):
    cdef c_StarWorld* _star_ptr   # non-owning; ownership via BaseWorld._world_ptr
    cpdef dict get_config_dict(self)
