# distutils: language = c++
"""
base.pxd
Cython declarations for TidalPy's base layer class (Phase 1).

Exports c_LayerEOSData, c_BaseLayerConfig, c_BaseLayer, and the Python
wrapper BaseLayer so other extensions can cimport and use C-speed access.

Usage::

    from TidalPy.structures_x.layers.base cimport (
        BaseLayer, c_BaseLayer, c_LayerEOSData)
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr

from TidalPy.Utilities_x.classes_x.classes cimport (
    TidalPyBaseClass,
    StructureBase,
    c_StructureBase,
)
from TidalPy.Material_x.eos.material_eos cimport c_MaterialEOSBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "eos_data_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayerEOSData:
        c_LayerEOSData()
        cpp_bool is_populated() const
        double   get_density(double radius_m) const
        double   get_gravity(double radius_m) const
        double   get_pressure(double radius_m) const
        void populate(
            const vector[double]& radius_m,
            const vector[double]& density_kgm3,
            const vector[double]& gravity_ms2,
            const vector[double]& pressure_pa)


cdef extern from "base_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_BaseLayerConfig:
        string   name
        int      layer_index
        double   radius_inner_m
        double   radius_outer_m
        double   mass_kg
        string   material_name
        cpp_bool is_tidal
        double   tidal_scale

    cdef cppclass c_BaseLayer(c_StructureBase):
        c_BaseLayer()
        c_BaseLayer(const c_BaseLayerConfig& config) except +
        const string& get_name()               const
        int      get_layer_index()             const
        double   get_radius_inner()            const
        double   get_radius_outer()            const
        double   get_thickness()               const
        double   get_volume()                  const
        double   get_surface_area_inner()      const
        double   get_surface_area_outer()      const
        const string& get_material_name()      const
        cpp_bool get_is_tidal()                const
        double   get_tidal_scale()             const
        cpp_bool get_eos_data_populated()      const
        double   get_density(double radius_m)  const
        double   get_gravity(double radius_m)  const
        double   get_pressure(double radius_m) const
        void     update_eos_data(const c_LayerEOSData& data)
        void     set_eos(unique_ptr[c_MaterialEOSBase] eos)
        c_MaterialEOSBase* get_eos() const
        cpp_bool get_eos_set() const
        cpp_bool get_viscoelastic_populated() const
        double   get_shear_modulus(double radius_m) const
        double   get_bulk_modulus(double radius_m) const
        double   get_shear_viscosity(double radius_m) const
        double   get_bulk_viscosity(double radius_m) const
        double   get_premelt_shear_modulus(double radius_m) const
        double   get_premelt_bulk_modulus(double radius_m) const
        double   get_premelt_shear_viscosity(double radius_m) const
        double   get_premelt_bulk_viscosity(double radius_m) const


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class BaseLayer(StructureBase):
    cdef unique_ptr[c_BaseLayer] _layer_ptr   # owns the most-derived C++ layer object
    cpdef dict get_config_dict(self)
