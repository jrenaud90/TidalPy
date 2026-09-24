# distutils: language = c++
"""Cython declarations for TidalPy's layered world class."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector
from libcpp.optional cimport optional

from CyRK cimport ODEMethod

from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport (
    BaseWorld, c_BaseWorld, c_WorldConfig, c_TideConfig, c_TideSolveConfig,
    c_Heating3DCollapseConfig, c_Heating3DCollapsed)
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage
from TidalPy.dynamics_x.spin cimport Spin, c_Spin
from TidalPy.Tides_x.love.love cimport c_LoveNumbers


cdef extern from "thermal_layout_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_LayerThermal:
        cpp_bool in_network
        cpp_bool boundary_fallback
        double temperature
        double top_temperature
        double node_temperature
        double heat_flow_in
        double heat_flow_out
        double boundary_thickness
        double rayleigh_number
        double nusselt_number
        double heating


cdef extern from "layered_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayerLove:
        size_t              layer_index
        double              tidal_scale
        c_LoveNumbers       love
        cpp_complex[double] shear_modulus
        double              volume

    cdef cppclass c_RadialSegment:
        size_t   world_layer
        double   radius_inner
        double   radius_outer
        cpp_bool molten

    cdef cppclass c_Grid3DAxes:
        const double* radii
        size_t        num_radii
        const double* colatitudes
        size_t        num_colatitudes
        const double* longitudes
        size_t        num_longitudes
        const double* times
        size_t        num_times

    cdef cppclass c_WorldEOSSolveConfig:
        double    surface_pressure
        size_t    slices_per_layer
        double    G_to_use
        ODEMethod integration_method
        double    rtol
        double    atol
        double    pressure_tol
        size_t    max_iters
        cpp_bool  nondimensionalize
        double    temperature
        cpp_bool  solve_temperature
        double    surface_temperature
        double    time
        size_t    max_thermal_passes
        double    thermal_tol
        double    radius_tol
        cpp_bool  reset_layer_masses
        cpp_bool  verbose

    cdef cppclass c_LoveSolveConfig:
        double    frequency
        int       degree_l
        vector[int] bc_models
        void      set_bc_models(const int* models_ptr, size_t num_models) except +
        int       love_method
        double    fixed_q
        double    fixed_dt
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

    cdef cppclass c_EOSSolverOverrides:
        optional[ODEMethod] integration_method
        optional[double]    rtol
        optional[double]    atol
        optional[double]    pressure_tol
        optional[size_t]    max_iters
        optional[size_t]    slices_per_layer
        optional[cpp_bool]  nondimensionalize
        optional[cpp_bool]  solve_temperature

    cdef cppclass c_RadialSolverOverrides:
        optional[ODEMethod] integration_method
        optional[double]    rtol
        optional[double]    atol
        optional[cpp_bool]  use_kamata
        optional[double]    start_radius_tol
        optional[cpp_bool]  scale_rtols
        optional[size_t]    max_num_steps
        optional[size_t]    expected_size
        optional[size_t]    max_ram_MB
        optional[cpp_bool]  nondimensionalize

    cdef cppclass c_LayeredWorld(c_BaseWorld):
        c_LayeredWorld()
        c_LayeredWorld(const c_WorldConfig& cfg) except +
        void         add_layer(unique_ptr[c_BaseLayer] layer) except +
        cpp_bool     accepts_layer(const c_BaseLayer& layer) except +
        string       layer_rejection_reason(const c_BaseLayer& layer) except +
        c_BaseLayer* get_layer(size_t index) except +
        void         update_after_layer_geometry_change() except +
        size_t       get_num_layers() const
        double       calc_total_mass() const
        double       calc_internal_heating(double time) const
        cpp_bool     validate_layers() const
        void         solve_eos(const c_WorldEOSSolveConfig& cfg) except +
        double       get_temperature(double radius)
        double       get_heat_flow(double radius)
        size_t       get_thermal_passes()
        cpp_bool     get_thermal_converged()
        cpp_bool     get_geometry_converged()
        double       calc_layer_temperature_rate(size_t layer_index)
        const vector[c_LayerThermal]& get_layer_thermal()
        double       get_density(double radius) const
        double       get_gravity(double radius) const
        double       get_pressure(double radius) const
        double       get_shear_modulus(double radius) const
        double       get_bulk_modulus(double radius) const
        double       get_shear_viscosity(double radius) const
        double       get_bulk_viscosity(double radius) const
        double       get_melt_fraction(double radius) const
        void         get_eos_state(double radius, double* y_out) const
        cpp_complex[double] calc_complex_shear_modulus(double radius, double frequency) const
        cpp_complex[double] calc_complex_bulk_modulus(double radius, double frequency) const
        cpp_bool     get_eos_solved() const
        cpp_bool     get_all_eos_set() const
        cpp_bool     get_eos_success() const
        const string& get_eos_message() const
        int          get_eos_iterations() const
        cpp_bool     get_eos_max_iters_hit() const
        double       get_eos_pressure_error() const
        double       get_surface_gravity_eos() const
        double       get_surface_pressure_eos() const
        double       get_central_pressure() const
        double       get_planet_mass_eos() const
        double       get_planet_moi_eos() const
        const vector[c_RadialSegment]& get_radial_segments() const
        vector[c_RadialSegment] get_molten_regions() const
        void         set_spin_model(const c_Spin& spin)
        const c_Spin& get_spin_model() const
        double       get_moment_of_inertia() const
        double       calc_spin_derivative(double host_mass) except +
        double       calc_synchronous_spin(double orbital_frequency) const
        const c_EOSSolution* get_eos_solution() const
        void                 set_eos_solver_overrides(const c_EOSSolverOverrides& overrides)
        void                 set_radial_solver_overrides(const c_RadialSolverOverrides& overrides)
        c_EOSSolverOverrides    get_eos_solver_overrides() const
        c_RadialSolverOverrides get_radial_solver_overrides() const
        c_WorldEOSSolveConfig make_eos_solve_config() const
        c_LoveSolveConfig    make_love_solve_config() const
        void                 solve_love_numbers(const c_LoveSolveConfig& cfg) except +
        void                 solve_love_numbers_supplied(
                const c_LoveSolveConfig& cfg,
                const cpp_complex[double]* shear_in,
                const cpp_complex[double]* bulk_in,
                const double* radius_in,
                size_t n_in) except +
        unique_ptr[c_RadialSolutionStorage] release_radial_storage() except +
        cpp_bool             get_love_solved() const
        cpp_bool             get_love_success() const
        int                  get_love_error_code() const
        const string&        get_love_message() const
        size_t               get_love_num_ytypes() const
        double               get_love_surface_amplification() const
        cpp_complex[double]  get_love_number_k(size_t ytype_idx) const
        cpp_complex[double]  get_love_number_h(size_t ytype_idx) const
        cpp_complex[double]  get_love_number_l(size_t ytype_idx) const
        cpp_complex[double]  get_love_surface_y(size_t ytype_idx, size_t y_idx) const
        cpp_complex[double]  get_radial_solution_y(double radius, size_t ytype_idx, size_t y_idx) const
        int                  get_love_method_last_int() const
        cpp_complex[double]  get_love_analytic_shear() const
        double               get_love_analytic_tidal_volume() const
        const vector[c_LayerLove]& get_love_layer_parts() const
        double               get_layer_tidal_scale(size_t index) except +
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
                                 double* out_heating,
                                 int num_threads) except +
        void                 get_3d_displacements_grid(
                                 const c_TideSolveConfig& state,
                                 const c_Grid3DAxes& axes,
                                 double* out_disp,
                                 int num_threads) except +
        void                 get_3d_stress_strain_grid(
                                 const c_TideSolveConfig& state,
                                 const c_Grid3DAxes& axes,
                                 double* out_stress,
                                 double* out_strain,
                                 int num_threads) except +
        c_Heating3DCollapsed calc_3d_tides(
                                 const c_TideSolveConfig& state,
                                 const double* radii,
                                 size_t num_radii,
                                 const double* colatitudes,
                                 size_t num_colatitudes,
                                 const double* longitudes,
                                 size_t num_longitudes,
                                 const double* times,
                                 size_t num_times,
                                 const c_Heating3DCollapseConfig& cfg) except +
        c_Heating3DCollapsed calc_3d_tides_layout(
                                 const double* radii,
                                 size_t num_radii,
                                 const double* colatitudes,
                                 size_t num_colatitudes,
                                 const double* longitudes,
                                 size_t num_longitudes,
                                 const double* times,
                                 size_t num_times,
                                 const c_Heating3DCollapseConfig& cfg) except +
        void                 calc_3d_tides_into(
                                 const c_TideSolveConfig& state,
                                 const double* radii,
                                 size_t num_radii,
                                 const double* colatitudes,
                                 size_t num_colatitudes,
                                 const double* longitudes,
                                 size_t num_longitudes,
                                 const double* times,
                                 size_t num_times,
                                 const c_Heating3DCollapseConfig& cfg,
                                 double* out_values,
                                 double* out_layer_totals) except +


cdef extern from "profile_world_.hpp" namespace "tidalpy" nogil:
    # Build a layered world whose layers interpolate their own slice of a radial profile. The one
    # implementation of that rule; `build_layered_world_from_profile` below and the standalone
    # `RadialSolver_x.radial_solver` both reach it, so neither can drift from the other.
    shared_ptr[c_LayeredWorld] c_build_world_from_layered_profile(
        const double* radius_ptr,
        const double* density_ptr,
        const double* shear_modulus_ptr,
        const double* bulk_modulus_ptr,
        size_t num_slices,
        const double* upper_radius_bylayer_ptr,
        const int* layer_type_ptr,
        const cpp_bool* is_static_ptr,
        const cpp_bool* is_incompressible_ptr,
        size_t num_layers,
        double planet_bulk_density,
        const string& name
    ) except +


cdef class LayeredWorld(BaseWorld):
    cdef c_LayeredWorld* _layered_ptr   # non-owning; ownership via BaseWorld._world_ptr
    @staticmethod
    cdef LayeredWorld _wrap(shared_ptr[c_BaseWorld] ptr)
    # Cached non-owning layer views, built once (lazily) and invalidated by add_layer so the
    # wrappers are not rebuilt on every world.<layer> / get_layer access.
    cdef list _layer_views
    cdef dict _layer_view_by_name
    # Weak references to every layer view handed out (cached or from add_layer), so a load that replaces the
    # layers can detach them instead of leaving them pointing at freed memory.
    cdef list _issued_views
    cdef void _track_view(self, BaseLayer view) except *
    # Scalar dispatch for the vectorized real-valued radius getters (nogil-callable
    # so the float-or-ndarray wrappers can loop without the GIL).
    cdef double _eval_real(self, int kind, double radius) noexcept nogil
    cpdef dict get_config_dict(self)
    cdef list _ensure_layer_views(self)
