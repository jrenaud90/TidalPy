from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex

from CyRK cimport ODEMethod
from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage

cdef extern from "solver_.hpp":
    int c_radial_solver(
        c_RadialSolutionStorage* solution_storage_ptr,
        size_t total_slices,
        double* radius_array_in_ptr,
        double* density_array_in_ptr,
        cpp_complex[double]* complex_bulk_modulus_in_ptr,
        cpp_complex[double]* complex_shear_modulus_in_ptr,
        double frequency,
        double planet_bulk_density,
        size_t num_layers,
        int* layer_types_ptr,
        cpp_bool* is_static_bylayer_ptr,
        cpp_bool* is_incompressible_bylayer_ptr,
        double surface_pressure,
        int degree_l,
        size_t num_bc_models,
        int* bc_models_ptr,
        int core_model,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        ODEMethod integration_method_int,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_bylayer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool nondimensionalize,
        cpp_bool use_prop_matrix,
        int* eos_integration_method_int_bylayer_ptr,
        ODEMethod eos_integration_method,
        double eos_rtol,
        double eos_atol,
        double eos_pressure_tol,
        int eos_max_iters,
        cpp_bool verbose,
        cpp_bool warnings
    ) noexcept nogil

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
