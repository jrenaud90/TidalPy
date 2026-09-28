# distutils: language = c++
"""Cython declarations for TidalPy's gas-giant world class."""

from libcpp.memory cimport shared_ptr

from TidalPy.Utilities.classes.classes cimport c_TidalPyBaseClass
from TidalPy.Structures.worlds.base cimport c_WorldConfig, c_BaseWorld
from TidalPy.Structures.worlds.layered cimport LayeredWorld, c_LayeredWorld


cdef extern from "gasgiant_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_GasGiantWorld(c_LayeredWorld):
        c_GasGiantWorld()
        c_GasGiantWorld(const c_WorldConfig& cfg) except +


cdef class GasGiantWorld(LayeredWorld):
    cdef c_GasGiantWorld* _gasgiant_ptr   # non-owning; ownership via BaseWorld._world_ptr
    cdef void _bind(self, shared_ptr[c_BaseWorld] ptr)
    @staticmethod
    cdef GasGiantWorld _wrap(shared_ptr[c_BaseWorld] ptr)
