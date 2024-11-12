from libcpp cimport bool as cpp_bool

from CyRK cimport CySolverResult

from TidalPy.RadialSolver.solutions cimport RadialSolutionStorageCC


cdef void cf_shooting_solver(
    RadialSolutionStorageCC* solution_storage_ptr,
    size_t total_slices,
    double* radius_array_ptr,
    double frequency,
    double planet_bulk_density,
    CySolverResult** eos_solution_bylayer_ptr,
    size_t num_layers,
    int* layer_types_ptr,
    int* is_static_by_layer_ptr,
    int* is_incompressible_by_layer_ptr,
    double* upper_radius_by_layer_ptr,
    size_t* first_slice_index_by_layer_ptr,
    size_t* num_slices_by_layer_ptr,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use = *,
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
    cpp_bool verbose = *,
    cpp_bool raise_on_fail = *
    ) noexcept nogil
