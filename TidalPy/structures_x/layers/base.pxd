# distutils: language = c++
"""Cython declarations for TidalPy's base layer: c_LayerEOSData, c_BaseLayerConfig, c_BaseLayer, and the Python
wrapper BaseLayer."""

from libc.stdint cimport uint32_t
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


cdef extern from "eos_data_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayerEOSData:
        c_LayerEOSData()
        cpp_bool is_populated() const
        double   get_density(double radius) const
        double   get_gravity(double radius) const
        double   get_pressure(double radius) const
        void populate(
            const vector[double]& radius,
            const vector[double]& density_kgm3,
            const vector[double]& gravity_ms2,
            const vector[double]& pressure) except +


cdef extern from "base_.hpp" namespace "tidalpy" nogil:
    const char* c_layer_class_name(uint32_t class_id)

    cdef cppclass c_BaseLayerConfig:
        string             name
        int                layer_index
        double             radius_inner
        double             radius_outer
        double             mass
        string             material_name
        cpp_bool           is_tidal
        cpp_bool           is_volume_fixed
        double             tidal_scale

    cdef cppclass c_BaseLayer(c_StructureBase):
        c_BaseLayer()
        c_BaseLayer(const c_BaseLayerConfig& config) except +
        const string& get_name()               const
        int      get_layer_index()             const
        double   get_radius_inner()            const
        double   get_radius_outer()            const
        double   get_thickness()               const
        double   get_volume()                  const
        double   get_density_bulk()            const
        double   get_surface_area_inner()      const
        double   get_surface_area_outer()      const
        const string& get_material_name()      const
        cpp_bool get_is_tidal()                const
        cpp_bool get_is_volume_fixed()         const
        void     set_is_volume_fixed(cpp_bool)
        void     set_radii(double radius_inner, double radius_outer)
        double   get_tidal_scale()             const
        void     set_tidal_scale(double tidal_scale)
        double   calc_tidal_scale(double planet_volume) const
        uint32_t get_layer_class_id()          const
        double   get_tidal_heating()           const
        cpp_bool get_eos_data_populated()      const
        double   get_density(double radius)  const
        double   get_gravity(double radius)  const
        double   get_pressure(double radius) const
        void     update_eos_data(const c_LayerEOSData& data)
        void     set_eos(unique_ptr[c_MaterialEOSBase] eos)
        c_MaterialEOSBase* get_eos() const
        cpp_bool get_eos_set() const
        cpp_bool get_viscoelastic_populated() const
        double   get_shear_modulus(double radius) const
        double   get_bulk_modulus(double radius) const
        double   get_shear_viscosity(double radius) const
        double   get_bulk_viscosity(double radius) const
        double   get_melt_fraction(double radius) const
        void     get_eos_state(double radius, double* y_out) const


cdef extern from "eos_layout_.hpp" nogil:
    const size_t C_EOS_DY_VALUES
    const size_t C_EOS_GRAVITY_INDEX
    const size_t C_EOS_PRESSURE_INDEX
    const size_t C_EOS_DENSITY_INDEX
    const size_t C_EOS_SHEAR_MODULUS_INDEX
    const size_t C_EOS_BULK_MODULUS_INDEX
    const size_t C_EOS_SHEAR_VISCOSITY_INDEX
    const size_t C_EOS_BULK_VISCOSITY_INDEX
    const size_t C_EOS_MELT_FRACTION_INDEX


# Fills the dense EOS layout at one radius for the object behind owner (a layer, a world).
ctypedef void (*cy_eos_state_fn)(const void* owner, double radius, double* y_out) noexcept nogil

cdef object cy_eos_fields(const void* owner, cy_eos_state_fn fill, object radius, tuple indices)


cdef class BaseLayer(StructureBase):
    cdef unique_ptr[c_BaseLayer] _layer_ptr   # owns the most-derived C++ layer object
    cdef cpp_bool _is_view                    # True => non-owning view into a world-owned layer
    cdef object   _world_ref                  # keep-alive ref to the owning world (views only)
    cpdef dict get_config_dict(self)
    # Scalar kernel behind the vectorized real-valued radius getters (see _apply_real in base.pyx).
    cdef double _eval_real(self, int kind, double radius) noexcept nogil
    # Initialize as a non-owning view; subclass `_view` factories set their own typed pointer first.
    cdef void _init_view(self, c_BaseLayer* ptr, object world)
    @staticmethod
    cdef BaseLayer _view(c_BaseLayer* ptr, object world)
