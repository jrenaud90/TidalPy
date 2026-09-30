# distutils: language = c++
"""Cython declarations for TidalPy's gas layer: c_GasConfig, c_GasLayer, and the Python wrapper GasLayer."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Structures.layers.base cimport BaseLayer, c_BaseLayer, c_BaseLayerConfig


cdef extern from "gas_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_GasConfig(c_BaseLayerConfig):
        double              mean_molecular_weight
        double              adiabatic_index
        double              reference_temperature
        double              reference_density

    cdef cppclass c_GasLayer(c_BaseLayer):
        c_GasLayer() except +
        c_GasLayer(const c_GasConfig& cfg) except +
        # Property getters
        double get_mean_molecular_weight()  const
        double get_adiabatic_index()        const
        double get_reference_temperature()  const
        double get_reference_density()      const


cdef class GasLayer(BaseLayer):
    cdef c_GasLayer* _gas_ptr   # non-owning; ownership via BaseLayer._layer_ptr
    
    @staticmethod
    cdef GasLayer _view(c_GasLayer* ptr, object world)
    cpdef dict get_config_dict(self)
