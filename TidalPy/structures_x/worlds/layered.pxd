# distutils: language = c++
"""
layered.pxd
Cython declarations for TidalPy's layered world class (Phase 8b).
"""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from CyRK cimport ODEMethod

from TidalPy.structures_x.worlds.base cimport (
    BaseWorld, c_BaseWorld, c_WorldConfig, c_TideConfig, c_TideSolveConfig)
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "layered_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_WorldEOSSolveConfig:
        double    surface_pressure
        size_t    slices_per_layer
        double    G_to_use
        ODEMethod integration_method
        double    rtol
        double    atol
        double    pressure_tol
        size_t    max_iters
        double    temperature
        cpp_bool  verbose

    cdef cppclass c_LoveSolveConfig:
        double    frequency_rad_s
        int       degree_l
        cpp_bool  solve_tidal
        cpp_bool  use_prop_matrix
        int       core_model
        cpp_bool  use_kamata
        cpp_bool  nondimensionalize
        double    starting_radius
        double    start_radius_tol
        ODEMethod integration_method
        double    rtol
        double    atol
        cpp_bool  scale_rtols
        size_t    max_num_steps
        size_t    expected_size
        size_t    max_ram_MB
        double    max_step
        cpp_bool  verbose
        cpp_bool  warnings
        double    eos_rtol
        double    eos_atol
        double    eos_pressure_tol
        int       eos_max_iters

    cdef cppclass c_LayeredWorld(c_BaseWorld):
        c_LayeredWorld()
        c_LayeredWorld(const c_WorldConfig& cfg) except +
        void         add_layer(unique_ptr[c_BaseLayer] layer) except +
        cpp_bool     accepts_layer(const c_BaseLayer& layer) const
        c_BaseLayer* get_layer(size_t index) except +
        size_t       get_num_layers() const
        double       calc_total_mass() const
        double       calc_internal_heating(double time_s) const
        cpp_bool     validate_layers() const
        void         solve_eos(const c_WorldEOSSolveConfig& cfg) except +
        double       get_density(double radius_m) const
        double       get_gravity(double radius_m) const
        double       get_pressure(double radius_m) const
        double       get_shear_modulus(double radius_m) const
        double       get_bulk_modulus(double radius_m) const
        double       get_shear_viscosity(double radius_m) const
        double       get_bulk_viscosity(double radius_m) const
        double       get_premelt_shear_modulus(double radius_m) const
        double       get_premelt_bulk_modulus(double radius_m) const
        double       get_premelt_shear_viscosity(double radius_m) const
        double       get_premelt_bulk_viscosity(double radius_m) const
        cpp_complex[double] calc_complex_shear_modulus(double radius_m, double frequency_rad_s) const
        cpp_complex[double] calc_complex_bulk_modulus(double radius_m, double frequency_rad_s) const
        cpp_bool     get_eos_solved() const
        cpp_bool     get_all_eos_set() const
        cpp_bool     get_eos_success() const
        const string& get_eos_message() const
        int          get_eos_iterations() const
        double       get_eos_pressure_error() const
        double       get_surface_gravity_eos() const
        double       get_surface_pressure_eos() const
        double       get_central_pressure() const
        double       get_planet_mass_eos() const
        double       get_planet_moi_eos() const
        const c_EOSSolution* get_eos_solution() const
        void                 solve_love_numbers(const c_LoveSolveConfig& cfg) except +
        void                 solve_love_numbers_supplied(const c_LoveSolveConfig& cfg, const cpp_complex[double]* shear_in, const cpp_complex[double]* bulk_in, const double* radius_in, size_t n_in) except +
        unique_ptr[c_RadialSolutionStorage] release_radial_storage()
        cpp_bool             get_love_solved() const
        cpp_bool             get_love_success() const
        int                  get_love_error_code() const
        const string&        get_love_message() const
        size_t               get_love_num_ytypes() const
        cpp_complex[double]  get_love_number_k(size_t ytype_idx) const
        cpp_complex[double]  get_love_number_h(size_t ytype_idx) const
        cpp_complex[double]  get_love_number_l(size_t ytype_idx) const
        cpp_complex[double]  get_love_surface_y(size_t ytype_idx, size_t y_idx) const
        # Global (1D) tidal dissipation: the model/config/result accessors are inherited from
        # c_BaseWorld; c_LayeredWorld only adds the rheology-capable calc_tides + layer heating.
        void                 calc_tides(const c_TideSolveConfig& state) except +
        double               get_layer_tidal_heating(size_t index) const
        # On-demand 3D tidal stress/strain/heating (rheology model; truncation from the tide config).
        double               get_3d_tidal_heating(
                                 const c_TideSolveConfig& state,
                                 double radius,
                                 double colatitude) except +
        void                 get_3d_tidal_heating_array(
                                 const c_TideSolveConfig& state,
                                 const double* radii,
                                 const double* colatitudes,
                                 size_t num_points,
                                 double* out_heating) except +


# =====================================================================================================================
# Cython wrapper class declaration
# =====================================================================================================================
cdef class LayeredWorld(BaseWorld):
    cdef c_LayeredWorld* _layered_ptr   # non-owning; ownership via BaseWorld._world_ptr
    # Cached non-owning layer views, built once (lazily) and invalidated by add_layer so the
    # wrappers are not rebuilt on every world.<layer> / get_layer access.
    cdef list _layer_views
    cdef dict _layer_view_by_name
    # Scalar dispatch for the vectorized real-valued radius getters (nogil-callable
    # so the float-or-ndarray wrappers can loop without the GIL).
    cdef double _eval_real(self, int kind, double radius_m) noexcept nogil
    cpdef dict get_config_dict(self)
    cdef list _ensure_layer_views(self)
