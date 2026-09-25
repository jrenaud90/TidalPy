# distutils: language = c++
"""Cython declarations for TidalPy's star world class."""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr, shared_ptr

from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig
from TidalPy.stellar_x.luminosity cimport LuminosityBase, c_LuminosityBase


cdef extern from "stellar_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_StarConfig(c_WorldConfig):
        double   effective_temperature
        double   luminosity

    cdef cppclass c_StarWorld(c_BaseWorld):
        c_StarWorld()
        c_StarWorld(const c_StarConfig& cfg) except +
        double get_effective_temperature() const
        double get_luminosity()            const
        double calc_luminosity_from_temperature(double temperature) const
        double calc_temperature_from_luminosity(double luminosity)  const
        void   set_effective_temperature(double temperature)
        void   set_luminosity(double luminosity)
        void   set_luminosity_model(unique_ptr[c_LuminosityBase] model)
        const c_LuminosityBase* get_luminosity_model() const
        cpp_bool has_luminosity_model() const
        double calc_luminosity_from_mass() except +
        double calc_effective_temperature_from_mass() except +
        void   update_luminosity_from_mass() except +


cdef class StarWorld(BaseWorld):
    cdef c_StarWorld* _star_ptr   # non-owning; ownership via BaseWorld._world_ptr
    cdef void _bind(self, shared_ptr[c_BaseWorld] ptr)
    cpdef dict get_config_dict(self)
    @staticmethod
    cdef StarWorld _wrap(shared_ptr[c_BaseWorld] ptr)
