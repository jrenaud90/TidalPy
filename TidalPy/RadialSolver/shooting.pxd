from libcpp cimport bool as cpp_bool
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

from TidalPy.RadialSolver.rs_solution cimport RadialSolutionStorageCC


cdef size_t cf_find_num_shooting_solutions(
    int layer_type,
    bint is_static,
    bint is_incompressible
    ) noexcept nogil


cdef int cf_shooting_solver(
    RadialSolutionStorageCC* solution_storage_ptr,
    double frequency,
    double planet_bulk_density,
    int* layer_types_ptr,
    bint* is_static_by_layer_ptr,
    bint* is_incompressible_by_layer_ptr,
    vector[size_t] first_slice_index_by_layer_vec,
    vector[size_t] num_slices_by_layer_vec,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use,
    int degree_l,
    cpp_bool use_kamata,
    double starting_radius,
    double start_radius_tolerance,
    int integration_method,
    double integration_rtol,
    double integration_atol,
    cpp_bool scale_rtols_by_layer_type,
    size_t max_num_steps,
    size_t expected_size,
    size_t max_ram_MB,
    double max_step,
    cpp_bool verbose
    ) noexcept nogil
