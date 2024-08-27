from libcpp cimport bool as cpp_bool

from TidalPy.RadialSolver.solutions cimport RadialSolutionStorageCC

cdef void cf_shooting_solver(
    RadialSolutionStorageCC* solution_storage_ptr,
    size_t total_slices,
    double* radius_array_ptr,
    double* density_array_ptr,
    double* gravity_array_ptr,
    double* bulk_modulus_array_ptr,
    double complex* complex_shear_modulus_array_ptr,
    double frequency,
    double planet_bulk_density,
    size_t num_layers,
    int* layer_types_ptr,
    int* is_static_by_layer_ptr,
    int* is_incompressible_by_layer_ptr,
    double* upper_radius_by_layer_ptr,
    size_t num_bc_models,
    int* bc_models_ptr,
    unsigned int degree_l = *,
    cpp_bool use_kamata = *,
    unsigned char integration_method = *,
    double integration_rtol = *,
    double integration_atol = *,
    cpp_bool scale_rtols_by_layer_type = *,
    size_t max_num_steps = *,
    size_t expected_size = *,
    size_t max_ram_MB = *,
    double max_step = *,
    cpp_bool limit_solution_to_radius = *,
    cpp_bool nondimensionalize = *,
    cpp_bool verbose = *,
    cpp_bool raise_on_fail = *
    ) noexcept
