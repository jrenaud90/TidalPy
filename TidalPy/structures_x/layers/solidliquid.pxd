# distutils: language = c++
"""Cython declarations for TidalPy's solid/liquid layer: c_SolidLiquidConfig, c_SolidLiquidLayer, and the Python
wrapper SolidLiquidLayer."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer, c_PhysicsConfig, c_BaseLayer
from TidalPy.cooling_x.cooling cimport c_CoolingBase
from TidalPy.radiogenics_x.radiogenics cimport c_RadiogenicsBase


cdef extern from "solidliquid_.hpp" namespace "tidalpy" nogil:

    # Adds no fields to c_PhysicsConfig.
    cdef cppclass c_SolidLiquidConfig(c_PhysicsConfig):
        pass

    cdef cppclass c_SolidLiquidLayer(c_PhysicsLayer):
        c_SolidLiquidLayer() except +
        c_SolidLiquidLayer(const c_SolidLiquidConfig& cfg) except +
        # Thermal constants of the material (read from the layer's EOS model)
        double get_thermal_conductivity() const
        double get_thermal_expansion()    const
        double get_heat_capacity()        const
        # Calculations
        double calc_thermal_conductivity(double temperature) const
        double calc_thermal_diffusivity(double temperature)  const
        double calc_adiabatic_temperature_gradient(double temperature, double pressure) const
        double calc_heat_flux_conductive(double temperature_base, double temperature_top) const
        double calc_radiogenic_heating(double time, double mass) const
        # Sub-model flags
        cpp_bool get_cooling_set()     const
        cpp_bool get_radiogenics_set() const
        c_CoolingBase*     get_cooling_model()     const
        c_RadiogenicsBase* get_radiogenics_model() const
        # Sub-model setters (transfer ownership)
        void     set_cooling(unique_ptr[c_CoolingBase] cooling)
        void     set_radiogenics(unique_ptr[c_RadiogenicsBase] radiogenics)


cdef class SolidLiquidLayer(PhysicsLayer):
    cdef c_SolidLiquidLayer* _solidliquid_ptr   # non-owning; ownership via BaseLayer._layer_ptr
    cpdef dict get_config_dict(self)
    
    @staticmethod
    cdef SolidLiquidLayer _view(c_SolidLiquidLayer* ptr, object world)
