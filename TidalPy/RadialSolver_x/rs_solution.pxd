cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr

from TidalPy.utilities.dimensions.nondimensional cimport NonDimensionalScalesCC
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution


cdef extern from "rs_solution_.hpp" nogil:

    cdef cppclass c_RadialSolutionStorage:

        c_RadialSolutionStorage()
        c_RadialSolutionStorage(
            size_t num_ytypes,
            double* upper_radius_bylayer_ptr,
            size_t num_layers,
            double* radius_array_ptr,
            size_t size_radius_array,
            int degree_l)

        cpp_bool success
        int error_code
        int degree_l
        cpp_string message
        size_t num_ytypes
        size_t num_slices
        size_t num_layers
        size_t total_size
        unique_ptr[c_EOSSolution] eos_solution_uptr
        vector[double] full_solution_vec
        vector[double] complex_love_vec
        vector[size_t] shooting_method_steps_taken_vec

        c_EOSSolution* get_eos_solution_ptr()
        void change_radius_array(
            double* new_radius_array_ptr,
            size_t new_size_radius_array,
            cpp_bool array_changed)
        void find_love()
        void dimensionalize_data(
            NonDimensionalScalesCC* nondim_scales,
            cpp_bool redimensionalize)


cdef class RadialSolverSolution:

    # Size and state information
    cdef size_t radius_array_size
    cdef public size_t num_ytypes
    cdef public size_t num_layers
    cdef cpp_bool ytype_names_set
    cdef char* ytypes[5]

    # Main storage container
    cdef unique_ptr[c_RadialSolutionStorage] solution_storage_uptr
    cdef c_RadialSolutionStorage* solution_storage_ptr

    # Result pointers and data
    cdef cnp.ndarray full_solution_arr

    # Love number information
    cdef cnp.ndarray complex_love_arr

    # EOS solution arrays
    cdef cnp.ndarray radius_array_cnp
    cdef cnp.ndarray gravity_array_cnp
    cdef cnp.ndarray pressure_array_cnp
    cdef cnp.ndarray mass_array_cnp
    cdef cnp.ndarray moi_array_cnp
    cdef cnp.ndarray density_array_cnp
    cdef cnp.ndarray shear_modulus_array_cnp
    cdef cnp.ndarray bulk_modulus_array_cnp

    # Shooting method diagnostics
    cdef cnp.ndarray shooting_method_steps_taken_array
    cdef cnp.ndarray eos_steps_taken_array

    cdef void finalize_python_storage(self) noexcept

    cdef void set_model_names(
        self,
        int* bc_models_ptr) noexcept nogil

    cdef void change_radius_array(
        self,
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        cpp_bool array_changed = *) noexcept
