# distutils: language = c++
"""
gasgiant.pxd
Cython declarations for TidalPy's gas-giant world class (Phase 8b).
"""

from TidalPy.structures_x.worlds.base cimport c_WorldConfig
from TidalPy.structures_x.worlds.layered cimport LayeredWorld, c_LayeredWorld


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "gasgiant_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_GasGiantWorld(c_LayeredWorld):
        c_GasGiantWorld()
        c_GasGiantWorld(const c_WorldConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class GasGiantWorld(LayeredWorld):
    cdef c_GasGiantWorld* _gasgiant_ptr   # non-owning; ownership via BaseWorld._world_ptr
