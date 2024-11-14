from libcpp cimport bool as cpp_bool
from libcpp.memory cimport shared_ptr

from TidalPy.RadialSolver.rs_solution cimport RadialSolutionStorageCC

cdef void cf_matrix_propagate(
    shared_ptr[RadialSolutionStorageCC] solution_storage_sptr,
    double frequency,
    double planet_bulk_density,
    # TODO: In the future the propagation matrix should take in layer types and multiple layers
    # int* layer_types_ptr,
    # int* is_static_by_layer_ptr,
    # int* is_incompressible_by_layer_ptr,
    size_t num_bc_models,
    int* bc_models_ptr,
    double G_to_use = *,
    unsigned int degree_l = *,
    unsigned char core_condition = *,
    cpp_bool verbose = *,
    cpp_bool raise_on_fail = *
    ) noexcept nogil
