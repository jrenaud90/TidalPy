
from libcpp cimport bool as cpp_bool

from CyRK.utils.memory cimport shared_ptr

from TidalPy.RadialSolver.rs_solution cimport RadialSolutionStorageCC


cdef void cf_radial_solver(
        shared_ptr[RadialSolutionStorageCC] solution_storage_sptr,
        double starting_radius,
        size_t total_slices,
        double* radius_array_in_ptr,
        double* density_array_in_ptr,
        double complex* complex_bulk_modulus_in_ptr,
        double complex* complex_shear_modulus_in_ptr,
        double frequency,
        double planet_bulk_density,
        size_t num_layers,
        int* layer_types_ptr,
        int* is_static_by_layer_ptr,
        int* is_incompressible_by_layer_ptr,
        double surface_pressure,
        unsigned int degree_l,
        size_t num_bc_models,
        int* bc_models_ptr,
        unsigned char core_condition,
        cpp_bool use_kamata,
        unsigned char integration_method_int,
        double integration_rtol,
        double integration_atol,
        cpp_bool scale_rtols_by_layer_type,
        size_t max_num_steps,
        size_t expected_size,
        size_t max_ram_MB,
        double max_step,
        cpp_bool nondimensionalize,
        cpp_bool use_prop_matrix,
        unsigned int eos_integration_method,
        double eos_rtol,
        double eos_atol,
        double eos_pressure_tol,
        unsigned int eos_max_iters,
        cpp_bool verbose,
        cpp_bool warnings,
        cpp_bool raise_on_fail,
        ) noexcept nogil
