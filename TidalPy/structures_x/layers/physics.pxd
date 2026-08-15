# distutils: language = c++
"""
physics.pxd
Cython declarations for TidalPy's physics layer class.

Exports c_PhysicsConfig, c_PhysicsLayer, and the Python wrapper PhysicsLayer
so other extensions can cimport and use C-speed access.

Usage::

    from TidalPy.structures_x.layers.physics cimport (
        PhysicsLayer, c_PhysicsLayer, c_PhysicsConfig)
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer, c_TidalScaleMethod
from TidalPy.Tides_x.love.love cimport c_LoveNumbers
from TidalPy.rheology_x.rheology cimport c_RheologyBase
from TidalPy.viscosity_x.viscosity cimport c_ViscosityBase
from TidalPy.partial_melt_x.partial_melt cimport c_PartialMeltBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "physics_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_PhysicsConfig:
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
        # PhysicsLayer additions:
        double              shear_modulus_static_pa
        double              bulk_modulus_static_pa
        double              shear_viscosity_static_pas
        double              bulk_viscosity_static_pas
        c_LoveNumbers       love_numbers
        # Radial-solver layer classification flags:
        cpp_bool            is_solid
        cpp_bool            is_static
        cpp_bool            is_incompressible

    cdef cppclass c_PhysicsLayer(c_BaseLayer):
        c_PhysicsLayer() except +
        c_PhysicsLayer(const c_PhysicsConfig& cfg) except +
        double              get_shear_modulus_static()               const
        double              get_bulk_modulus_static()                const
        double              get_shear_viscosity_static()             const
        double              get_bulk_viscosity_static()              const
        c_LoveNumbers       get_love_numbers()                       const
        cpp_complex[double] get_love_number_k()                      const
        cpp_complex[double] get_love_number_h()                      const
        cpp_complex[double] get_love_number_l()                      const
        double              calc_tidal_susceptibility()              const
        cpp_complex[double] calc_complex_shear_modulus(double freq)  const
        cpp_complex[double] calc_complex_bulk_modulus(double freq)   const
        cpp_bool            get_shear_rheology_set()                 const
        cpp_bool            get_bulk_rheology_set()                  const
        void                set_shear_rheology(unique_ptr[c_RheologyBase] shear)
        void                set_bulk_rheology(unique_ptr[c_RheologyBase] bulk)
        void                set_shear_viscosity(unique_ptr[c_ViscosityBase] viscosity)
        void                set_bulk_viscosity(unique_ptr[c_ViscosityBase] viscosity)
        void                set_partial_melt(unique_ptr[c_PartialMeltBase] partial_melt)
        cpp_bool            get_shear_viscosity_set()                const
        cpp_bool            get_bulk_viscosity_set()                 const
        cpp_bool            get_partial_melt_set()                   const
        cpp_bool            get_is_solid()                           const
        cpp_bool            get_is_static()                          const
        cpp_bool            get_is_incompressible()                  const
        void                set_is_solid(cpp_bool)
        void                set_is_static(cpp_bool)
        void                set_is_incompressible(cpp_bool)


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class PhysicsLayer(BaseLayer):
    cdef c_PhysicsLayer* _physics_ptr   # non-owning; ownership via BaseLayer._layer_ptr
    cpdef dict get_config_dict(self)
    
    @staticmethod
    cdef PhysicsLayer _view(c_PhysicsLayer* ptr, object world)
