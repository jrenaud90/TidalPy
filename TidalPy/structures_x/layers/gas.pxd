# distutils: language = c++
"""
gas.pxd
Cython declarations for TidalPy's gas layer class.

Exports c_GasConfig, c_GasLayer, and the Python wrapper GasLayer so other
extensions can cimport and use C-speed access.

Usage::

    from TidalPy.structures_x.layers.gas cimport (
        GasLayer, c_GasLayer, c_GasConfig)
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.complex cimport complex as cpp_complex

from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer, c_BaseLayer
from TidalPy.structures_x.layers.base cimport c_TidalScaleMethod
from TidalPy.Tides_x.love.love cimport c_LoveNumbers


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "gas_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_GasConfig:
        # Inherited from c_BaseLayerConfig:
        string              name
        int                 layer_index
        double              radius_inner
        double              radius_outer
        double              mass
        string              material_name
        cpp_bool            is_tidal
        double              tidal_scale
        c_TidalScaleMethod  tidal_scale_method
        # From c_PhysicsConfig:
        double              shear_modulus_static
        double              bulk_modulus_static
        double              shear_viscosity_static
        double              bulk_viscosity_static
        c_LoveNumbers       love_numbers
        # GasLayer additions:
        double              mean_molecular_weight
        double              adiabatic_index
        double              reference_temperature
        double              reference_density

    cdef cppclass c_GasLayer(c_PhysicsLayer):
        c_GasLayer() except +
        c_GasLayer(const c_GasConfig& cfg) except +
        # Property getters
        double get_mean_molecular_weight()  const
        double get_adiabatic_index()        const
        double get_reference_temperature()  const
        double get_reference_density()      const
        # Calculations
        double calc_adiabatic_lapse_rate(double gravity)                      const
        double calc_scale_height(double temperature, double gravity)        const
        double calc_pressure_ideal_gas(double temperature, double density) const
        double calc_sound_speed(double temperature)                               const


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class GasLayer(PhysicsLayer):
    cdef c_GasLayer* _gas_ptr   # non-owning; ownership via BaseLayer._layer_ptr
    
    @staticmethod
    cdef GasLayer _view(c_GasLayer* ptr, object world)
    cpdef dict get_config_dict(self)
