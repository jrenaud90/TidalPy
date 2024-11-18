from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr

cimport numpy as cnp
cnp.import_array()

from CyRK cimport CySolverResult

cdef extern from "eos_solution_.cpp" nogil:

    const size_t EOS_Y_VALUES
    const size_t EOS_EXTRA_VALUES
    const size_t EOS_DY_VALUES

    cdef cppclass EOSSolutionCC:
            int error_code
            int iterations
            cpp_bool success
            cpp_bool max_iters_hit
            cpp_bool radius_array_set
            cpp_bool other_vecs_set

            char* message_ptr

            size_t current_layers_saved
            size_t num_layers
            size_t radius_array_size

            double pressure_error
            double surface_gravity
            double surface_pressure
            double central_pressure
            double radius
            double mass
            double moi

            vector[double] upper_radius_bylayer_vec
            vector[shared_ptr[CySolverResult]] cysolver_results_sptr_bylayer_vec

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
                const size_t num_layers,
                double* radius_array_ptr,
                const size_t radius_array_size)
            
            void save_cyresult(
                shared_ptr[CySolverResult] new_cysolver_result_sptr)
            
            void call(
                const size_t layer_index,
                const double radius,
                double* y_interp_ptr)
            void call_vectorize(
                const size_t layer_index,
                const double* radius_array_ptr,
                size_t len_radius_array,
                double* y_interp_ptr)

            void change_radius_array(
                double* radius_array_ptr,
                const size_t radius_array_size)
            
            void interpolate_full_planet()
