# distutils: language = c++
"""
layered.pxd
Cython declarations for TidalPy's layered world class (Phase 8b).
"""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr

from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "layered_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayeredWorld(c_BaseWorld):
        c_LayeredWorld()
        c_LayeredWorld(const c_WorldConfig& cfg) except +
        void         add_layer(unique_ptr[c_BaseLayer] layer) except +
        cpp_bool     accepts_layer(const c_BaseLayer& layer) const
        c_BaseLayer* get_layer(size_t index) except +
        size_t       get_num_layers() const
        double       calc_total_mass() const
        double       calc_internal_heating(double time_s) const
        cpp_bool     validate_layers() const


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class LayeredWorld(BaseWorld):
    cdef c_LayeredWorld* _layered_ptr   # non-owning; ownership via BaseWorld._world_ptr
    cpdef dict get_config_dict(self)
