from libcpp cimport bool as cpp_bool

from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage


cdef extern from "matrix_.hpp" nogil:
    int c_matrix_propagate(
        c_RadialSolutionStorage* solution_storage_ptr,
        double frequency,
        double planet_bulk_density,
        size_t* first_slice_index_by_layer_ptr,
        size_t* num_slices_by_layer_ptr,
        size_t num_layers_for_slices,
        size_t num_bc_models,
        int* bc_models_ptr,
        double G_to_use,
        int degree_l,
        double starting_radius,
        double start_radius_tolerance,
        int core_model,
        cpp_bool verbose
        ) noexcept
