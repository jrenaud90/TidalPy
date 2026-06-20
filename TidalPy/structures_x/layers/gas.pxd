# distutils: language = c++
"""
gas.pxd
Cython declarations for TidalPy's gas layer class (Phase 4).

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
        double              radius_inner_m
        double              radius_outer_m
        double              mass_kg
        string              material_name
        cpp_bool            is_tidal
        double              tidal_scale
        c_TidalScaleMethod  tidal_scale_method
        # From c_PhysicsConfig:
        double              shear_modulus_static_pa
        double              bulk_modulus_static_pa
        double              shear_viscosity_static_pas
        double              bulk_viscosity_static_pas
        c_LoveNumbers       love_numbers
        # GasLayer additions:
        double              mean_molecular_weight_kg_mol
        double              adiabatic_index
        double              reference_temperature_k
        double              reference_density_kg_m3

    cdef cppclass c_GasLayer(c_PhysicsLayer):
        c_GasLayer() except +
        c_GasLayer(const c_GasConfig& cfg) except +
        # Property getters
        double get_mean_molecular_weight()  const
        double get_adiabatic_index()        const
        double get_reference_temperature()  const
        double get_reference_density()      const
        # Calculations
        double calc_adiabatic_lapse_rate(double gravity_m_s2)                      const
        double calc_scale_height(double temperature_k, double gravity_m_s2)        const
        double calc_pressure_ideal_gas(double temperature_k, double density_kg_m3) const
        double calc_sound_speed(double temperature_k)                               const


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class GasLayer(PhysicsLayer):
    cdef c_GasLayer* _gas_ptr   # non-owning; ownership via BaseLayer._layer_ptr
    cpdef dict get_config_dict(self)
