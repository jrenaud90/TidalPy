from libcpp cimport bool as cpp_bool
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex

from CyRK cimport ODEMethod
from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage


# The world types are redeclared here rather than cimported from structures_x.worlds.layered, which declares
# the same C++ classes. A cimport would be compile-time for these plain C++ declarations, but Cython also emits
# a runtime import of every extension type the cimported .pxd declares (LayeredWorld, and through its own
# cimports BaseWorld and BaseLayer). RadialSolver_x/__init__ imports this module and layered.pyx cimports
# RadialSolver_x.rs_solution, so that runtime import would close an import cycle. Only the members this module
# calls are declared; the header itself is the single definition, so a signature that changes there fails to
# compile here rather than drifting silently.
cdef extern from "layered_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_WorldEOSSolveConfig:
        double    surface_pressure
        size_t    slices_per_layer
        ODEMethod integration_method
        double    rtol
        double    atol
        double    pressure_tol
        size_t    max_iters
        cpp_bool  nondimensionalize

    cdef cppclass c_LoveSolveConfig:
        double    frequency
        int       degree_l
        void      set_bc_models(const int* models_ptr, size_t num_models) except +
        int       love_method
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

    cdef cppclass c_LayeredWorld:
        c_WorldEOSSolveConfig make_eos_solve_config() const
        void solve_eos(const c_WorldEOSSolveConfig& cfg) except +
        cpp_bool get_eos_max_iters_hit() const
        void solve_love_numbers_supplied(
            const c_LoveSolveConfig& cfg,
            const cpp_complex[double]* shear_in,
            const cpp_complex[double]* bulk_in,
            const double* radius_in,
            size_t n_in) except +
        unique_ptr[c_RadialSolutionStorage] release_radial_storage()


cdef extern from "profile_world_.hpp" namespace "tidalpy":
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
        const cpp_string& name
    ) except +

cdef extern from "solver_.hpp":
    void c_validate_and_prep_radial_inputs(
        size_t total_slices,
        const double* radius_array,
        const double* density_array,
        double frequency,
        size_t num_layers,
        const vector[cpp_string]& layer_types,
        const cpp_bool* is_static_bylayer,
        const cpp_bool* is_incompressible_bylayer,
        const double* upper_radius_bylayer_array,
        cpp_bool use_prop_matrix,
        double starting_radius,
        const vector[cpp_string]& solve_for,
        const cpp_string& integration_method,
        const vector[cpp_string]& eos_method_bylayer,
        const cpp_string& eos_integration_method,
        cpp_bool warnings,
        int* layer_types_out_ptr,
        int* bc_models_out_ptr,
        size_t& num_bc_models_out,
        ODEMethod& integration_method_out,
        vector[int]& eos_integration_method_int_bylayer_out,
        ODEMethod& eos_integration_method_out
    ) except +
