from libcpp cimport bool as cpp_bool
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

from TidalPy.RadialSolver.rs_solution cimport RadialSolutionStorageCC

cdef int cf_matrix_propagate(
    RadialSolutionStorageCC* solution_storage_ptr,
    double frequency,
    double planet_bulk_density,
    # TODO: In the future the propagation matrix should take in layer types and multiple layers
    # int* layer_types_ptr,
    # int* is_static_by_layer_ptr,
    # int* is_incompressible_by_layer_ptr,
    vector[size_t] first_slice_index_by_layer_vec,
    vector[size_t] num_slices_by_layer_vec,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use,
    int degree_l,
    double starting_radius,
    double start_radius_tolerance,
    int core_model,
    cpp_bool verbose
    ) noexcept nogil
