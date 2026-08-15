# distutils: language = c++
"""
stellar.pxd
Cython declarations for TidalPy's star world class.
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr

from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld
from TidalPy.stellar_x.luminosity cimport LuminosityBase, c_LuminosityBase


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
        void   set_luminosity_model(unique_ptr[c_LuminosityBase] model)
        cpp_bool has_luminosity_model() const
        double calc_luminosity_from_mass() except +
        double calc_effective_temperature_from_mass() except +
        void   update_luminosity_from_mass() except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class StarWorld(BaseWorld):
    cdef c_StarWorld* _star_ptr   # non-owning; ownership via BaseWorld._world_ptr
    cpdef dict get_config_dict(self)
    @staticmethod
    cdef StarWorld _wrap(shared_ptr[c_BaseWorld] ptr)
