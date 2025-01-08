
from libcpp cimport bool as cpp_bool
from libcpp.memory cimport shared_ptr

from TidalPy.RadialSolver.rs_solution cimport RadialSolutionStorageCC


cdef int cf_radial_solver(
    shared_ptr[RadialSolutionStorageCC] solution_storage_sptr,
    size_t total_slices,
    double* radius_array_in_ptr,
    double* density_array_in_ptr,
    double complex* complex_bulk_modulus_in_ptr,
    double complex* complex_shear_modulus_in_ptr,
    double frequency,
    double planet_bulk_density,
    size_t num_layers,
    int* layer_types_ptr,
    bint* is_static_bylayer_ptr,
    bint* is_incompressible_bylayer_ptr,
    double surface_pressure,
    int degree_l,
    size_t num_bc_models,
    int* bc_models_ptr,
    int core_model,
    cpp_bool use_kamata,
    double starting_radius,
    double start_radius_tolerance,
    int integration_method_int,
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
    int eos_integration_method,
    double eos_rtol,
    double eos_atol,
    double eos_pressure_tol,
    int eos_max_iters,
    cpp_bool verbose,
    cpp_bool warnings
    ) noexcept nogil
