# distutils: language = c++
"""Cython declarations for TidalPy's gas-giant world class."""

from libcpp.memory cimport shared_ptr

from TidalPy.Structures.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig


cdef extern from "gasgiant_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_GasGiantWorld(c_BaseWorld):
        c_GasGiantWorld()
        c_GasGiantWorld(const c_WorldConfig& cfg) except +


cdef class GasGiantWorld(BaseWorld):
    @staticmethod
    cdef GasGiantWorld _wrap(shared_ptr[c_BaseWorld] ptr)
