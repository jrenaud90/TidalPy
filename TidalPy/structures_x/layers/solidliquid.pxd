# distutils: language = c++
"""
solidliquid.pxd
Cython declarations for TidalPy's solid/liquid layer class (Phase 3).

Exports c_SolidLiquidConfig, c_SolidLiquidLayer, and the Python wrapper
SolidLiquidLayer so other extensions can cimport and use C-speed access.

Usage::

    from TidalPy.structures_x.layers.solidliquid cimport (
        SolidLiquidLayer, c_SolidLiquidLayer, c_SolidLiquidConfig)
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer, c_BaseLayer
from TidalPy.structures_x.layers.base cimport c_TidalScaleMethod
from TidalPy.Tides_x.love.love cimport c_LoveNumbers
from TidalPy.cooling_x.cooling cimport c_CoolingBase
from TidalPy.radiogenics_x.radiogenics cimport c_RadiogenicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "solidliquid_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_SolidLiquidConfig:
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
        # SolidLiquidLayer additions:
        double              thermal_conductivity_ref_w_mk
        double              thermal_expansion_ref_1_k
        double              heat_capacity_ref_j_kgk
        double              activation_energy_j_mol
        double              activation_volume_m3_mol
        double              solidus_temperature_k
        double              liquidus_temperature_k
        double              melt_fraction_exponent
        double              reference_density_kg_m3
        double              reference_temperature_k
        double              melt_viscosity_reduction

    cdef cppclass c_SolidLiquidLayer(c_PhysicsLayer):
        c_SolidLiquidLayer() except +
        c_SolidLiquidLayer(const c_SolidLiquidConfig& cfg) except +
        # Thermal property getters
        double get_thermal_conductivity_ref()   const
        double get_thermal_expansion_ref()      const
        double get_heat_capacity_ref()          const
        double get_activation_energy()          const
        double get_activation_volume()          const
        double get_solidus_temperature()        const
        double get_liquidus_temperature()       const
        double get_melt_fraction_exponent()     const
        double get_reference_density()          const
        double get_reference_temperature()      const
        double get_melt_viscosity_reduction()   const
        # Calculations
        double calc_melt_fraction(double temperature_k, double pressure_pa)     const
        double calc_viscosity(double temperature_k, double pressure_pa)         const
        double calc_shear_modulus(double temperature_k, double pressure_pa)     const
        double calc_thermal_conductivity(double temperature_k)                  const
        double calc_thermal_diffusivity(double temperature_k)                   const
        double calc_adiabatic_temperature_gradient(double temperature_k, double pressure_pa) const
        double calc_heat_flux_conductive(double temperature_base_k, double temperature_top_k) const
        double calc_radiogenic_heating(double time_s, double mass_kg)           const
        # Sub-model flags
        cpp_bool get_cooling_set()      const
        cpp_bool get_radiogenics_set()  const
        # Sub-model setters (transfer ownership)
        void     set_cooling(unique_ptr[c_CoolingBase] cooling)
        void     set_radiogenics(unique_ptr[c_RadiogenicsBase] radiogenics)


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class SolidLiquidLayer(PhysicsLayer):
    cdef c_SolidLiquidLayer* _solidliquid_ptr   # non-owning; ownership via BaseLayer._layer_ptr
    cpdef dict get_config_dict(self)
    
    @staticmethod
    cdef SolidLiquidLayer _view(c_SolidLiquidLayer* ptr, object world)
