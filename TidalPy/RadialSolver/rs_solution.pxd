cimport numpy as np

from libcpp cimport bool as cpp_bool

from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr

from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC

# Need to include love_.cpp and eos_solution_.cpp in order to get solutions.cpp to see it and use it to link
cdef extern from "love_.cpp" nogil:
    pass

cdef extern from "eos_solution_.cpp" nogil:
    pass

cdef extern from "solutions_.cpp" nogil:
    
    const int MAX_NUM_Y
    const int MAX_NUM_Y_REAL

    cdef cppclass RadialSolutionStorageCC:
        
        cpp_bool success
        char num_ytypes
        char* message_ptr
        size_t num_slices
        size_t num_layers
        size_t total_size

        RadialSolutionStorageCC()
        RadialSolutionStorageCC(
            size_t num_slices,
            char num_ytypes
            )
        size_t num_slices
        size_t total_size
        char[256] message
        char* message_ptr
        cpp_bool success
        char num_ytypes

        shared_ptr[EOSSolutionCC] eos_solution_sptr
        vector[double] full_solution_vec
        vector[double] complex_love_vec

        double* full_solution_ptr
        double* complex_love_ptr
        double* radius_array_ptr
        double* gravity_array_ptr
        double* pressure_array_ptr
        double* mass_array_ptr
        double* moi_array_ptr
        double* density_array_ptr
        double* complex_shear_array_ptr
        double* complex_bulk_array_ptr

        void set_message(const char* new_message)
        void find_love()


cdef class RadialSolverSolution:

    # Size and state information
    cdef size_t num_slices
    cdef char num_ytypes
    cdef cpp_bool ytype_names_set
    cdef char* ytypes[5]

    # Main storage container
    cdef shared_ptr[RadialSolutionStorageCC] solution_storage_sptr

    # Result pointers and data
    cdef np.ndarray full_solution_arr

    # Love number information
    cdef np.ndarray complex_love_arr

    # EOS solution arrays
    cdef np.ndarray gravity_array
    cdef np.ndarray pressure_array
    cdef np.ndarray density_array
    cdef np.ndarray shear_modulus_array
    cdef np.ndarray bulk_modulus_array

    cdef void set_model_names(self, int* bc_models_ptr) noexcept nogil


cdef size_t cf_find_num_shooting_solutions(
    int layer_type,
    bint is_static,
    bint is_incompressible
    ) noexcept nogil
