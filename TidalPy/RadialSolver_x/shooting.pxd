from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector

from CyRK cimport ODEMethod

from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage

cdef extern from "shooting_.hpp" nogil:
    inline size_t c_find_num_shooting_solutions(
        int layer_type,
        cpp_bool is_static,
        cpp_bool is_incompressible
    ) noexcept

    int c_shooting_solver(
        c_RadialSolutionStorage* solution_storage_ptr,
        double frequency,
        double planet_bulk_density,
        int* layer_types_ptr,
        cpp_bool* is_static_by_layer_ptr,
        cpp_bool* is_incompressible_by_layer_ptr,
        vector[size_t]& first_slice_index_by_layer_vec,
        vector[size_t]& num_slices_by_layer_vec,
        size_t num_bc_models,
        int* bc_models_ptr,
        double G_to_use,
        int degree_l,
        cpp_bool use_kamata,
        double starting_radius,
        double start_radius_tolerance,
        ODEMethod integration_method,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_by_layer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool verbose
    ) noexcept