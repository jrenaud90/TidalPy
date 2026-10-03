# distutils: language = c++
"""Cython declarations for TidalPy's base world class: c_WorldConfig, c_BaseWorld, and the Python wrapper
BaseWorld."""

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector
from libcpp.optional cimport optional

from CyRK cimport ODEMethod

from TidalPy.Utilities.classes.classes cimport (
    TidalPyBaseClass,
    StructureBase,
    c_StructureBase,
    c_TidalPyBaseClass,
    c_PhysicsBase,
)
from TidalPy.Tides.classes.tide cimport c_TideBase
from TidalPy.Structures.layers.layer cimport Layer, c_Layer
from TidalPy.Material.eos.eos_solution cimport c_EOSSolution, c_EOSZone
from TidalPy.RadialSolver.rs_solution cimport c_RadialSolutionStorage
from TidalPy.Dynamics.spin cimport Spin, c_Spin
from TidalPy.Tides.love.love cimport c_LoveNumbers


# =====================================================================================================================
# Global (1D) and 3D tidal solve config structs (global namespace, tide_result_.hpp)
# =====================================================================================================================
cdef extern from "tide_result_.hpp" nogil:
    cdef cppclass c_TideConfig:
        int min_degree_l
        int max_degree_l
        int eccentricity_truncation
        double eccentricity_exact_tolerance
        int obliquity_truncation
        cpp_bool layer_tidal_heating
        int love_method
        double love_fixed_q
        double love_fixed_dt

    cdef cppclass c_TideSolveConfig:
        double orbital_frequency
        double spin_frequency
        double eccentricity
        double obliquity
        double semi_major_axis
        double host_mass

    cdef cppclass c_Heating3DCollapseConfig:
        cpp_bool orbit_averaged
        cpp_bool latitude_summed
        cpp_bool longitude_summed
        cpp_bool radial_summed
        int latitude_nodes
        int longitude_nodes
        int radial_slices
        int num_threads
        cpp_bool latitude_analytic
        double colatitude_min
        double colatitude_max

    cdef cppclass c_Heating3DCollapsed:
        vector[double] values
        vector[size_t] shape
        vector[double] radii
        vector[double] colatitudes
        vector[double] longitudes
        vector[double] times
        vector[double] layer_totals
        size_t n_layers
        size_t n_times
        cpp_bool all_spatial_summed


cdef extern from "thermal_layout_.hpp" namespace "tidalpy" nogil:
    cdef cppclass c_LayerThermal:
        cpp_bool in_network
        cpp_bool boundary_fallback
        double temperature
        double top_temperature
        double base_temperature
        double node_temperature
        double heat_flow_in
        double heat_flow_out
        double boundary_thickness
        double rayleigh_number
        double nusselt_number
        cpp_bool magma_ocean
        cpp_bool boundary_fallback
        double reference_pressure
        double reference_viscosity
        double reference_melt_fraction
        double heating
        # One entry per heat source, in c_HeatSourceKind order (radiogenic, tidal, prescribed).
        double heating_by_source[3]
        double latent_capacity
        double thermal_capacity


cdef extern from "base_.hpp" namespace "tidalpy" nogil:
    cdef enum class c_HeatSourceKind:
        Radiogenic
        Tidal
        Prescribed

    cdef cppclass c_PrescribedLayerHeating:
        double power
        double specific_rate

    cdef cppclass c_WorldConfig:
        string name
        string world_type_str
        double radius
        double mass
        double albedo
        double emissivity
        double obliquity
        double spin_frequency

    cdef cppclass c_WorldEOSReport:
        cpp_bool solved
        cpp_bool success
        string message
        int iterations
        cpp_bool max_iters_hit
        double pressure_error
        double surface_gravity
        double surface_pressure
        double central_pressure
        double planet_mass
        double planet_moi
        size_t thermal_passes
        cpp_bool thermal_converged
        vector[double] radius
        vector[double] gravity
        vector[double] pressure
        vector[double] mass
        vector[double] moi
        vector[double] density
        vector[double] temperature
        vector[double] heat_flow
        vector[c_LayerThermal] layer_thermal
        vector[double] layer_temperature_rate
        vector[double] layer_radius_outer
        vector[c_EOSZone] zones

    cdef cppclass c_LayerLove:
        size_t layer_index
        double tidal_scale
        c_LoveNumbers love
        cpp_complex[double] shear_modulus
        double volume

    cdef cppclass c_Grid3DAxes:
        const double* radii
        size_t num_radii
        const double* colatitudes
        size_t num_colatitudes
        const double* longitudes
        size_t num_longitudes
        const double* times
        size_t num_times

    cdef cppclass c_WorldEOSSolveConfig:
        double surface_pressure
        size_t slices_per_layer
        double G_to_use
        ODEMethod integration_method
        double rtol
        double atol
        double pressure_tol
        size_t max_iters
        cpp_bool nondimensionalize
        double temperature
        cpp_bool solve_temperature
        double surface_temperature
        double time
        size_t max_thermal_passes
        double thermal_tol
        cpp_bool reset_layer_masses
        cpp_bool verbose

    cdef cppclass c_LoveSolveConfig:
        double frequency
        int degree_l
        vector[int] bc_models
        void set_bc_models(const int* models_ptr, size_t num_models) except +
        int love_method
        double fixed_q
        double fixed_dt
        int core_model
        cpp_bool use_kamata
        cpp_bool nondimensionalize
        double starting_radius
        double start_radius_tol
        ODEMethod integration_method
        double rtol
        double atol
        cpp_bool scale_rtols
        size_t max_num_steps
        size_t expected_size
        size_t max_ram_MB
        double max_step
        cpp_bool verbose
        cpp_bool warnings

    # One section's pinned solver keys (c_SolverOverrides), by the TidalPy_Configs.toml key; a value is a double
    # whatever its kind (kind: 0 an ODE method's enum value, 1 a real number, 2 a count, 3 a switch).
    cdef cppclass c_EOSSolverOverrides:
        void set(const string& key, double value) except +
        cpp_bool has(const string& key) except +
        double get(const string& key) except +
        @staticmethod
        vector[string] keys()
        @staticmethod
        int kind(const string& key) except +

    cdef cppclass c_RadialSolverOverrides:
        void set(const string& key, double value) except +
        cpp_bool has(const string& key) except +
        double get(const string& key) except +
        @staticmethod
        vector[string] keys()
        @staticmethod
        int kind(const string& key) except +

    cdef cppclass c_BaseWorld(c_StructureBase):
        c_BaseWorld()
        c_BaseWorld(const c_WorldConfig& cfg) except +

        # Identity, orbital and thermal scalars, and bulk geometry
        const string& get_name() const
        const string& get_world_type() const
        double get_albedo() const
        double get_emissivity() const
        double get_obliquity() const
        double get_spin_frequency() const
        double calc_surface_gravity() const
        double calc_escape_velocity() const
        double calc_mean_density() const
        double calc_equilibrium_temperature(double insolation_flux) const
        void set_name(const string& name) except +
        void set_spin_frequency(double freq)
        void set_obliquity(double obliq)

        # Layers
        void add_layer(unique_ptr[c_Layer] layer) except +
        string layer_rejection_reason(const c_Layer& layer) except +
        c_Layer* get_layer(size_t index) except +
        void update_after_layer_change() except +
        size_t get_num_layers() const
        double calc_total_mass() const
        double calc_internal_heating(double time) const
        cpp_bool validate_layers() const

        # EOS solve and profile reads
        void solve_eos(const c_WorldEOSSolveConfig& cfg) except +
        c_WorldEOSReport solve_eos_report(const c_WorldEOSSolveConfig& cfg) except +
        c_WorldEOSReport get_eos_report() except +
        double get_temperature(double radius)
        double get_heat_flow(double radius)
        size_t get_thermal_passes()
        cpp_bool get_thermal_converged()
        double calc_layer_temperature_rate(size_t layer_index)
        const vector[c_LayerThermal]& get_layer_thermal()
        double get_density(double radius) const
        double get_gravity(double radius) const
        double get_pressure(double radius) const
        double get_shear_modulus(double radius) const
        double get_bulk_modulus(double radius) const
        double get_shear_viscosity(double radius) const
        double get_bulk_viscosity(double radius) const
        double get_melt_fraction(double radius) const
        void get_eos_state(double radius, double* y_out) const
        # Vectorized profile reads; each takes the world's call lock once for the whole array.
        void get_eos_fields(
            const size_t* field_indices,
            size_t num_fields,
            const double* radii,
            size_t num_radii,
            double* values_out) const
        cpp_complex[double] calc_complex_shear_modulus(double radius, double frequency) const
        cpp_complex[double] calc_complex_bulk_modulus(double radius, double frequency) const
        void calc_complex_moduli(
            cpp_bool is_shear,
            const double* radii,
            size_t num_radii,
            double frequency,
            cpp_complex[double]* moduli_out) const
        cpp_bool get_eos_solved() const
        cpp_bool get_all_materials_set() const
        cpp_bool get_eos_success() const
        string get_eos_message() except +
        int get_eos_iterations() const
        cpp_bool get_eos_max_iters_hit() const
        double get_eos_pressure_error() const
        double get_surface_gravity_eos() const
        double get_surface_pressure_eos() const
        double get_central_pressure() const
        double get_planet_mass_eos() const
        double get_planet_moi_eos() const
        vector[c_EOSZone] get_zones_copy() except +
        vector[c_EOSZone] get_molten_regions() except +
        vector[double] get_tidal_heat_source() except +
        void clear_tidal_heating()
        void set_prescribed_heating(size_t layer_index, double power, double specific_rate) except +
        c_PrescribedLayerHeating get_prescribed_heating(size_t layer_index)
        double get_heating(double radius) noexcept
        const c_EOSSolution* get_eos_solution() const

        # Spin
        void set_spin_model(const c_Spin& spin)
        const c_Spin& get_spin_model() const
        double get_moment_of_inertia() const
        double get_moment_of_inertia_factor() const
        double calc_spin_derivative(double host_mass) except +
        double calc_synchronous_spin(double orbital_frequency) const

        # Pinned solver settings
        void set_eos_solver_overrides(const c_EOSSolverOverrides& overrides)
        void set_radial_solver_overrides(const c_RadialSolverOverrides& overrides)
        c_EOSSolverOverrides get_eos_solver_overrides() const
        c_RadialSolverOverrides get_radial_solver_overrides() const
        c_WorldEOSSolveConfig make_eos_solve_config() const
        c_LoveSolveConfig make_love_solve_config() const

        # Love-number solve
        void solve_love_numbers(const c_LoveSolveConfig& cfg) except +
        void solve_love_numbers_supplied(
            const c_LoveSolveConfig& cfg,
            const cpp_complex[double]* shear_in,
            const cpp_complex[double]* bulk_in,
            const double* radius_in,
            size_t n_in) except +
        unique_ptr[c_RadialSolutionStorage] release_radial_storage() except +
        cpp_bool get_love_solved() const
        cpp_bool get_love_success() const
        int get_love_error_code() const
        string get_love_message() except +
        size_t get_love_num_ytypes() const
        double get_love_surface_amplification() const
        double get_love_surface_rcond() const
        cpp_complex[double] get_love_number_k(size_t ytype_idx) const
        cpp_complex[double] get_love_number_h(size_t ytype_idx) const
        cpp_complex[double] get_love_number_l(size_t ytype_idx) const
        double get_love_q_k(size_t ytype_idx) const
        double get_love_lag_k(size_t ytype_idx) const
        cpp_complex[double] get_love_surface_y(size_t ytype_idx, size_t y_idx) const
        cpp_complex[double] get_radial_solution_y(double radius, size_t ytype_idx, size_t y_idx) const
        int get_love_method_last_int() const
        cpp_complex[double] get_love_analytic_shear() const
        double get_love_analytic_tidal_volume() const
        vector[c_LayerLove] get_love_layer_parts() except +

        # Global (1D) tidal dissipation; calc_tides is defined in world_tides_.hpp.
        void set_tide_model_handle(const shared_ptr[c_PhysicsBase]& model) except +
        cpp_bool get_tide_model_set() const
        const c_TideBase* get_tide_model() const
        shared_ptr[c_PhysicsBase] share_tide_model() except +
        void set_tide_config(const c_TideConfig& cfg) except +
        const c_TideConfig& get_tide_config() const
        void calc_tides(const c_TideSolveConfig& state) except +
        cpp_bool get_tide_state(c_TideSolveConfig& state_out) except +
        cpp_bool get_tides_solved() const
        double get_tidal_heating() const
        double get_tidal_heat_flux() const
        double get_tidal_dU_dM() const
        double get_tidal_dU_dw() const
        double get_tidal_dU_dO() const
        double get_tidal_dU_dM_minus_dw() const
        int get_num_tidal_modes() const
        cpp_complex[double] get_tidal_love_k(int degree_l, int m, int p, int q) const
        double get_layer_tidal_scale(size_t index) except +
        double get_layer_tidal_heating(size_t index) const

        # On-demand 3D tidal stress, strain, and heating (rheology model; truncation from the tide config).
        double get_3d_tidal_heating(
            const c_TideSolveConfig& state,
            double radius,
            double colatitude) except +
        void get_3d_tidal_heating_array(
            const c_TideSolveConfig& state,
            const double* radii,
            const double* colatitudes,
            size_t num_points,
            double* out_heating,
            int num_threads) except +
        void get_3d_displacements_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_disp,
            int num_threads) except +
        void get_3d_stress_strain_grid(
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
        void calc_3d_tides_into(
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


cdef extern from "factory_.hpp" namespace "tidalpy" nogil:
    # A world of its own concrete class rebuilt from one complete world record held in memory (a copy or a pickle).
    shared_ptr[c_BaseWorld] c_world_from_binary_bytes(const string& record_bytes, cpp_bool force) except +


cdef extern from "profile_world_.hpp" namespace "tidalpy" nogil:
    # Build a world whose layers interpolate their own slice of a radial profile. The one implementation of that
    # rule; `build_layered_world_from_profile` and the standalone `RadialSolver.radial_solver` both reach it, so
    # neither can drift from the other.
    shared_ptr[c_BaseWorld] c_build_world_from_layered_profile(
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


# The orbital state of one tide solve, from the arguments of calc_tides in their order.
cdef inline c_TideSolveConfig cy_tide_state(
        double orbital_frequency,
        double spin_frequency,
        double eccentricity,
        double obliquity,
        double semi_major_axis,
        double host_mass) noexcept nogil:
    cdef c_TideSolveConfig state
    state.orbital_frequency = orbital_frequency
    state.spin_frequency    = spin_frequency
    state.eccentricity      = eccentricity
    state.obliquity         = obliquity
    state.semi_major_axis   = semi_major_axis
    state.host_mass         = host_mass
    return state


# Fill the fields every world config shares from the constructor arguments. A property given as None takes the
# [worlds] default of the TidalPy configuration for the world type ([worlds.<type>] winning), as the world builder
# does, and the C++ default where the configuration has none. Returns those defaults, for the properties a constructor
# sets after building the world (cy_set_default_spin).
cdef inline dict cy_fill_world_config(
        c_WorldConfig* config,
        str name,
        double radius,
        double mass,
        str world_type,
        object albedo,
        object emissivity,
        object obliquity,
        object spin_frequency):
    # Deferred: the configs package imports the world classes.
    from TidalPy.Structures.configs.toml_loader import world_type_defaults
    cdef dict defaults = world_type_defaults(world_type)
    config.name           = name.encode("utf-8")
    config.world_type_str = world_type.encode("utf-8")
    config.radius         = radius
    config.mass           = mass
    if albedo is None:
        albedo = defaults.get("albedo", config.albedo)
    if emissivity is None:
        emissivity = defaults.get("emissivity", config.emissivity)
    if obliquity is None:
        obliquity = defaults.get("obliquity_rad", config.obliquity)
    if spin_frequency is None:
        spin_frequency = defaults.get("spin_frequency_rad_s", config.spin_frequency)
    config.albedo         = <double>albedo
    config.emissivity     = <double>emissivity
    config.obliquity      = <double>obliquity
    config.spin_frequency = <double>spin_frequency
    return defaults


cdef class BaseWorld(StructureBase):
    # Owns the most-derived C++ world object (shared so a System can co-own it).
    cdef shared_ptr[c_BaseWorld] _world_ptr
    # The normalized config the world was built from, or None.
    cdef public dict source_config
    # A data-file world's config as given, for save_to_toml, or None.
    cdef public dict portable_config
    # The world's get_config_dict() at the end of its build, which save_to_toml compares the live state against, or
    # None.
    cdef public dict built_config
    # Cached non-owning layer views, built once (lazily) and invalidated by add_layer so the wrappers are not rebuilt
    # on every world.<layer> / get_layer access.
    cdef list _layer_views
    cdef dict _layer_view_by_name
    # Weak references to every layer view handed out (cached or from add_layer), so a load that replaces the layers
    # can detach them instead of leaving them pointing at freed memory.
    cdef list _issued_views
    cpdef dict get_config_dict(self)
    # Point this wrapper at a C++ world; each subclass also sets its own typed pointer.
    cdef void _bind(self, shared_ptr[c_BaseWorld] ptr)
    # Wrap an already-constructed C++ world without building a new one; each subclass returns its own type.
    @staticmethod
    cdef BaseWorld _wrap(shared_ptr[c_BaseWorld] ptr)
    cdef void _track_view(self, Layer view) except *
    cdef list _ensure_layer_views(self)


# The Python-side configurations of a world (source, portable, built) that a copy carries over, and their setter.
cdef dict cy_world_configs(BaseWorld world)
cdef void cy_set_world_configs(BaseWorld world, dict configs) except *


# A newly built world's spin model takes the [worlds] moment-of-inertia factor (cy_fill_world_config's defaults).
cdef inline void cy_set_default_spin(BaseWorld world, dict defaults) except *:
    if "moment_of_inertia_factor" in defaults:
        world.set_spin_model(Spin(moment_of_inertia_factor=defaults["moment_of_inertia_factor"]))
