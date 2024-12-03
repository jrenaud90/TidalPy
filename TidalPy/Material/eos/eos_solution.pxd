from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr

cimport numpy as cnp
cnp.import_array()

from CyRK cimport CySolverResult

from TidalPy.utilities.dimensions.nondimensional cimport NonDimensionalScalesCC

cdef extern from "nondimensional_.hpp" nogil:
    pass

cdef extern from "eos_solution_.cpp" nogil:

    const size_t EOS_Y_VALUES
    const size_t EOS_EXTRA_VALUES
    const size_t EOS_DY_VALUES

    cdef cppclass EOSSolutionCC:
            int error_code
            int iterations
            int nondim_status
            int solution_nondim_status
            cpp_bool success
            cpp_bool max_iters_hit
            cpp_bool radius_array_set
            cpp_bool other_vecs_set

            char* message_ptr

            size_t current_layers_saved
            size_t num_layers
            size_t radius_array_size
            size_t num_cyolver_calls

            double pressure_error
            double surface_gravity
            double surface_pressure
            double central_pressure
            double radius
            double mass
            double moi

            double redim_length_scale
            double redim_gravity_scale
            double redim_mass_scale
            double redim_density_scale
            double redim_moi_scale
            double redim_pascal_scale

            vector[double] upper_radius_bylayer_vec
            vector[shared_ptr[CySolverResult]] cysolver_results_sptr_bylayer_vec
            vector[size_t] steps_taken_vec

            vector[double] radius_array_vec
            vector[double] gravity_array_vec
            vector[double] pressure_array_vec
            vector[double] mass_array_vec
            vector[double] moi_array_vec
            vector[double] density_array_vec

            # These complex arrays are stored as double arrays with twice the length (Cython and C++ don't play nicely with complex across all systems)
            vector[double] complex_shear_array_vec
            vector[double] complex_bulk_array_vec

            EOSSolutionCC()
            EOSSolutionCC(
                    double* upper_radius_bylayer_ptr,
                    size_t num_layers,
                    double* radius_array_ptr,
                    size_t radius_array_size
                )
            
            void save_cyresult(
                shared_ptr[CySolverResult] new_cysolver_result_sptr)
            void save_steps_taken(
                size_t steps_taken)
            
            void call(
                const size_t layer_index,
                const double radius,
                double* y_interp_ptr)

            void change_radius_array(
                double* new_radius_ptr,
                size_t new_radius_size)
            
            void interpolate_full_planet()

            void dimensionalize_data(
                NonDimensionalScalesCC* nondim_scales,
                cpp_bool redimensionalize)
