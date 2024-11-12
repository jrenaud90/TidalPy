from libcpp cimport bool as cpp_bool

from TidalPy.Material.eos.solver cimport EOSSolutionVec
from TidalPy.RadialSolver.solutions cimport RadialSolutionStorageCC

cdef void cf_matrix_propagate(
    RadialSolutionStorageCC* solution_storage_ptr,
    size_t total_slices,
    double* radius_array_ptr,
    double frequency,
    double planet_bulk_density,
    EOSSolutionVec* eos_solution_bylayer_ptr,
    size_t num_layers,
    # TODO: In the future the propagation matrix should take in layer types and multiple layers
    # int* layer_types_ptr,
    # int* is_static_by_layer_ptr,
    # int* is_incompressible_by_layer_ptr,
    double* upper_radius_by_layer_ptr,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use = *,
    unsigned int degree_l = *,
    unsigned char core_condition = *,
    cpp_bool verbose = *,
    cpp_bool raise_on_fail = *
    ) noexcept nogil
