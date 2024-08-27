from libcpp cimport bool as cpp_bool

from TidalPy.RadialSolver.solutions cimport RadialSolutionStorageCC

cdef void cf_matrix_propagate(
    RadialSolutionStorageCC* solution_storage_ptr,
    size_t total_slices,
    double* radius_array_ptr,
    double* density_array_ptr,
    double* gravity_array_ptr,
    double* bulk_modulus_ptr,
    double complex* complex_shear_modulus_array_ptr,
    double frequency,
    double planet_bulk_density,
    # TODO: In the future the propagation matrix should take in layer types and multiple layers
    # size_t num_layers,
    # int* layer_types_ptr,
    # int* is_static_by_layer_ptr,
    # int* is_incompressible_by_layer_ptr,
    # double* upper_radius_by_layer_ptr,
    size_t num_bc_models,
    int* bc_models_ptr,
    unsigned int degree_l = *,
    unsigned char core_condition = *,
    cpp_bool nondimensionalize = *,
    cpp_bool verbose = *,
    cpp_bool raise_on_fail = *
    ) noexcept nogil