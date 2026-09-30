# distutils: language = c++
"""Cython declarations for TidalPy's terrestrial world class."""

from libcpp.memory cimport shared_ptr

from TidalPy.Structures.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig


cdef extern from "terrestrial_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_TerrestrialWorld(c_BaseWorld):
        c_TerrestrialWorld()
        c_TerrestrialWorld(const c_WorldConfig& cfg) except +


cdef class TerrestrialWorld(BaseWorld):
    @staticmethod
    cdef TerrestrialWorld _wrap(shared_ptr[c_BaseWorld] ptr)
