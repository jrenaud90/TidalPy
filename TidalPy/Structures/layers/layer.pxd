# distutils: language = c++
"""Cython declarations for TidalPy's layer: c_LayerEOSData, c_LayerConfig, c_Layer, and the Python wrapper Layer."""

from libc.stdint cimport uint8_t
from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities.classes.classes cimport (
    StructureBase,
    c_StructureBase,
    c_PhysicsBase,
    c_ThermoPoint,
)
from TidalPy.Material.material cimport c_MaterialSwitches, c_MaterialState
from TidalPy.Cooling.cooling cimport c_CoolingBase
from TidalPy.Radiogenics.radiogenics cimport c_RadiogenicsBase


cdef extern from "eos_data_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayerEOSData:
        c_LayerEOSData()
        cpp_bool is_populated() const
        void populate(
            const vector[double]& radius,
            const vector[double]& density_kgm3,
            const vector[double]& gravity_ms2,
            const vector[double]& pressure) except +


cdef extern from "layer_.hpp" namespace "tidalpy" nogil:

    cdef enum class c_LayerState(uint8_t):
        Auto
        Solid
        Liquid

    const char* c_layer_state_name(c_LayerState state)
    c_LayerState c_layer_state_from_name(const string& name) except +

    cdef cppclass c_LayerConfig:
        string             name
        int                layer_index
        double             radius_inner
        double             radius_outer
        double             mass
        cpp_bool           use_tides
        cpp_bool           is_volume_fixed
        double             tidal_scale
        c_LayerState       state
        cpp_bool           is_static
        cpp_bool           is_incompressible
        double             temperature
        c_MaterialSwitches switches
        cpp_bool           use_heating

    cdef cppclass c_Layer(c_StructureBase):
        c_Layer()
        c_Layer(const c_LayerConfig& config) except +
        # Geometry and identification
        const string& get_name() const
        int      get_layer_index() const
        double   get_radius_inner() const
        double   get_radius_outer() const
        double   get_radius_mid() const
        double   get_thickness() const
        double   get_volume() const
        double   get_density_bulk() const
        double   get_surface_area_inner() const
        double   get_surface_area_outer() const
        cpp_bool get_is_volume_fixed() const
        void     set_is_volume_fixed(cpp_bool) except +
        void     set_radii(double radius_inner, double radius_outer) except +
        # Tides
        cpp_bool get_use_tides() const
        void     set_use_tides(cpp_bool)
        double   get_tidal_scale() const
        void     set_tidal_scale(double tidal_scale)
        double   calc_tidal_scale(double planet_volume) const
        double   get_tidal_heating() const
        # Material, switches, and state
        cpp_bool get_material_set() const
        void     set_material_model(const shared_ptr[c_PhysicsBase]& model) except +
        shared_ptr[c_PhysicsBase] share_material_model() const
        const c_MaterialSwitches& get_switches() const
        void     set_switches(const c_MaterialSwitches& switches) except +
        double   get_temperature() const
        void     set_temperature(double) except +
        cpp_bool get_use_heating() const
        void     set_use_heating(cpp_bool) except +
        void     calc_state(const c_ThermoPoint& point, c_MaterialState& out) const
        c_LayerState get_state() const
        void     set_state(c_LayerState state) except +
        cpp_bool get_is_liquid() const
        cpp_bool get_can_change_state() const
        cpp_bool get_is_static() const
        cpp_bool get_is_incompressible() const
        void     set_is_static(cpp_bool) except +
        void     set_is_incompressible(cpp_bool) except +
        # Rheology
        void     set_shear_rheology_model(const shared_ptr[c_PhysicsBase]& model) except +
        void     set_bulk_rheology_model(const shared_ptr[c_PhysicsBase]& model) except +
        shared_ptr[c_PhysicsBase] share_shear_rheology_model(cpp_bool is_override) const
        shared_ptr[c_PhysicsBase] share_bulk_rheology_model(cpp_bool is_override) const
        cpp_complex[double] calc_complex_shear_modulus(double frequency) const
        cpp_complex[double] calc_complex_bulk_modulus(double frequency) const
        cpp_complex[double] calc_complex_shear_modulus(double radius, double frequency) const
        cpp_complex[double] calc_complex_bulk_modulus(double radius, double frequency) const
        # Vectorized radius-resolved form; takes the owning world's call lock once for the whole array.
        void calc_complex_moduli(
            cpp_bool is_shear,
            const double* radii,
            size_t num_radii,
            double frequency,
            cpp_complex[double]* moduli_out) const
        # Cooling and radiogenics (shared)
        void set_cooling_model(const shared_ptr[c_PhysicsBase]& model) except +
        shared_ptr[c_PhysicsBase] share_cooling_model() const
        void set_radiogenics_model(const shared_ptr[c_PhysicsBase]& model) except +
        shared_ptr[c_PhysicsBase] share_radiogenics_model() const
        const c_CoolingBase*     get_cooling_model() const
        const c_RadiogenicsBase* get_radiogenics_model() const
        double calc_radiogenic_heating(double time, double mass) const
        # Solved profile
        cpp_bool get_eos_data_populated() const
        void     update_eos_data(const c_LayerEOSData& data) except +
        cpp_bool get_viscoelastic_populated() const
        # Vectorized profile read; takes the owning world's call lock once for the whole array.
        void get_eos_fields(
            const size_t* field_indices,
            size_t num_fields,
            const double* radii,
            size_t num_radii,
            double* values_out) const


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


cdef class Layer(StructureBase):
    cdef unique_ptr[c_Layer] _layer_ptr   # owns the C++ layer object
    cdef cpp_bool _is_view                # True => non-owning view into a world-owned layer
    cdef object   _world_ref              # keep-alive ref to the owning world (views only)
    cdef cpp_bool _detached               # True => a world load replaced the layer this view pointed at
    cdef object   __weakref__             # lets the owning world track the views it hands out
    cdef void _check_ptr(self) except *
    # Called by the owning world when a load replaces its layers: forget the C++ layer without deleting it.
    cdef void _detach(self) noexcept
    # Tell the owning world a view moved its layer's radii.
    cdef void _notify_world_of_move(self) except *
    cpdef dict get_config_dict(self)
    # Initialize as a non-owning view.
    cdef void _init_view(self, c_Layer* ptr, object world)
    @staticmethod
    cdef Layer _view(c_Layer* ptr, object world)
