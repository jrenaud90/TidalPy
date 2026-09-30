# distutils: language = c++
"""Cython declarations for TidalPy's base layer: c_LayerEOSData, c_BaseLayerConfig, c_BaseLayer, and the Python
wrapper BaseLayer."""

from libc.stdint cimport uint32_t
from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities.classes.classes cimport (
    TidalPyBaseClass,
    StructureBase,
    c_StructureBase,
)
from TidalPy.Material.eos.material_eos cimport c_MaterialEOSBase
from TidalPy.Rheology.rheology cimport c_RheologyBase
from TidalPy.Viscosity.viscosity cimport c_ViscosityBase
from TidalPy.PartialMelt.partial_melt cimport c_PartialMeltBase


cdef extern from "eos_data_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayerEOSData:
        c_LayerEOSData()
        cpp_bool is_populated() const
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
        # Radial-solver layer classification flags
        cpp_bool           is_solid
        cpp_bool           is_static
        cpp_bool           is_incompressible
        double             temperature
        cpp_bool           use_thermal_eos
        cpp_bool           use_heating

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
        void     set_is_volume_fixed(cpp_bool) except +
        void     set_radii(double radius_inner, double radius_outer)
        double   get_tidal_scale()             const
        void     set_tidal_scale(double tidal_scale)
        double   calc_tidal_scale(double planet_volume) const
        uint32_t get_layer_class_id()          const
        double   get_tidal_heating()           const
        cpp_bool get_eos_data_populated()      const
        void     update_eos_data(const c_LayerEOSData& data)
        void     set_eos(unique_ptr[c_MaterialEOSBase] eos) except +
        c_MaterialEOSBase* get_eos() const
        cpp_bool get_eos_set() const
        cpp_bool get_viscoelastic_populated() const
        double   get_shear_modulus(double radius) const
        void     get_eos_state(double radius, double* y_out) const
        # Vectorized profile read; takes the owning world's call lock once for the whole array.
        void     get_eos_fields(
                     const size_t* field_indices,
                     size_t num_fields,
                     const double* radii,
                     size_t num_radii,
                     double* values_out) const
        double              get_shear_modulus_static()               const
        double              get_bulk_modulus_static()                const
        double              get_shear_viscosity_static()             const
        double              get_bulk_viscosity_static()              const
        cpp_complex[double] calc_complex_shear_modulus(double freq)  const
        cpp_complex[double] calc_complex_bulk_modulus(double freq)   const
        cpp_complex[double] calc_complex_shear_modulus(double radius, double freq) const
        cpp_complex[double] calc_complex_bulk_modulus(double radius, double freq)  const
        # Vectorized radius-resolved form; takes the owning world's call lock once for the whole array.
        void                calc_complex_moduli(
                                cpp_bool is_shear,
                                const double* radii,
                                size_t num_radii,
                                double frequency,
                                cpp_complex[double]* moduli_out) const
        cpp_bool            get_shear_rheology_set()                 const
        c_RheologyBase*     get_shear_rheology_model()               const
        c_RheologyBase*     get_bulk_rheology_model()                const
        c_ViscosityBase*    get_shear_viscosity_model()              const
        c_ViscosityBase*    get_bulk_viscosity_model()               const
        c_PartialMeltBase*  get_partial_melt_model()                 const
        cpp_bool            get_bulk_rheology_set()                  const
        void                set_shear_rheology(unique_ptr[c_RheologyBase] shear) except +
        void                set_bulk_rheology(unique_ptr[c_RheologyBase] bulk) except +
        void                set_shear_viscosity(unique_ptr[c_ViscosityBase] viscosity) except +
        void                set_bulk_viscosity(unique_ptr[c_ViscosityBase] viscosity) except +
        void                set_partial_melt(unique_ptr[c_PartialMeltBase] partial_melt) except +
        cpp_bool            get_shear_viscosity_set()                const
        cpp_bool            get_bulk_viscosity_set()                 const
        cpp_bool            get_partial_melt_set()                   const
        cpp_bool            get_is_solid()                           const
        cpp_bool            get_is_static()                          const
        cpp_bool            get_is_incompressible()                  const
        void                set_is_solid(cpp_bool) except +
        void                set_is_static(cpp_bool) except +
        void                set_is_incompressible(cpp_bool) except +
        double              get_temperature()                          const
        cpp_bool            get_use_thermal_eos()                      const
        void                set_temperature(double) except +
        void                set_use_thermal_eos(cpp_bool) except +
        cpp_bool            get_use_heating()                          const
        void                set_use_heating(cpp_bool) except +



cdef extern from "eos_layout_.hpp" nogil:
    const size_t C_EOS_DY_VALUES
    const size_t C_EOS_GRAVITY_INDEX
    const size_t C_EOS_PRESSURE_INDEX
    const size_t C_EOS_DENSITY_INDEX
    const size_t C_EOS_SHEAR_MODULUS_INDEX
    const size_t C_EOS_BULK_MODULUS_INDEX
    const size_t C_EOS_SHEAR_VISCOSITY_INDEX
    const size_t C_EOS_BULK_VISCOSITY_INDEX
    const size_t C_EOS_TEMPERATURE_INDEX
    const size_t C_EOS_HEAT_FLOW_INDEX
    const size_t C_EOS_MELT_FRACTION_INDEX


# Fills entries field_indices[0 .. num_fields) of the dense EOS layout at radii[0 .. num_radii) for the object behind
# owner (a layer, a world), field-major: values_out[field_i * num_radii + radius_i]. The C++ call holds the object's
# call lock for the whole array, so it runs without the GIL and never waits for it while holding the lock.
ctypedef void (*cy_eos_fields_fn)(
    const void* owner,
    const size_t* field_indices,
    size_t num_fields,
    const double* radii,
    size_t num_radii,
    double* values_out) noexcept nogil

cdef object cy_eos_field(const void* owner, cy_eos_fields_fn fill, object radius, size_t field_index)
cdef object cy_eos_fields(const void* owner, cy_eos_fields_fn fill, object radius, tuple indices)

# Fills the c_BaseLayerConfig fields shared by every layer constructor; a layer subclass passes its own config.
cdef int cy_fill_base_layer_config(
    c_BaseLayerConfig* config,
    str name,
    int layer_index,
    double radius_inner,
    double radius_outer,
    double mass,
    str material_name,
    cpp_bool is_tidal,
    cpp_bool is_volume_fixed,
    object tidal_scale,
    cpp_bool is_solid,
    cpp_bool is_static,
    cpp_bool is_incompressible,
    double temperature,
    cpp_bool use_thermal_eos,
    cpp_bool use_heating) except -1


cdef class BaseLayer(StructureBase):
    cdef unique_ptr[c_BaseLayer] _layer_ptr   # owns the most-derived C++ layer object
    cdef cpp_bool _is_view                    # True => non-owning view into a world-owned layer
    cdef object   _world_ref                  # keep-alive ref to the owning world (views only)
    cdef cpp_bool _detached                   # True => a world load replaced the layer this view pointed at
    cdef object   __weakref__                 # lets the owning world track the views it hands out
    cdef void _check_ptr(self) except *
    # Called by the owning world when a load replaces its layers: forget the C++ layer without deleting it.
    cdef void _detach(self) noexcept
    # Tell the owning world a view moved its layer's radii.
    cdef void _notify_world_of_move(self) except *
    cpdef dict get_config_dict(self)
    # Initialize as a non-owning view; subclass `_view` factories set their own typed pointer first.
    cdef void _init_view(self, c_BaseLayer* ptr, object world)
    @staticmethod
    cdef BaseLayer _view(c_BaseLayer* ptr, object world)
