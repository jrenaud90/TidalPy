cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool

from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, unique_ptr

from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC

# Need to include love_.cpp and eos_solution_.cpp in order to get solutions.cpp to see it and use it to link
cdef extern from "love_.cpp" nogil:
    pass

cdef extern from "eos_solution_.cpp" nogil:
    pass

cdef extern from "rs_solution_.cpp" nogil:
    
    const int MAX_NUM_Y
    const int MAX_NUM_Y_REAL

    cdef cppclass RadialSolutionStorageCC:
        
        RadialSolutionStorageCC()
        RadialSolutionStorageCC(
            char num_ytypes,
            double* upper_radius_bylayer_ptr,
            size_t num_layers,
            double* radius_array_ptr,
            size_t size_radius_array
            )

        cpp_bool success
        int error_code
        char* message_ptr
        size_t num_ytypes
        size_t num_slices
        size_t num_layers
        size_t total_size
        unique_ptr[EOSSolutionCC] eos_solution_uptr
        vector[double] full_solution_vec
        vector[double] complex_love_vec

        EOSSolutionCC* get_eos_solution_ptr()
        void change_radius_array(
            double* new_radius_array_ptr,
            size_t new_size_radius_array,
            cpp_bool array_changed)
        void set_message(const char* new_message)
        void find_love()


cdef class RadialSolverSolution:

    # Size and state information
    cdef size_t radius_array_size
    cdef size_t num_ytypes
    cdef cpp_bool ytype_names_set
    cdef char* ytypes[5]

    # Main storage container
    cdef shared_ptr[RadialSolutionStorageCC] solution_storage_sptr

    # Result pointers and data
    cdef cnp.ndarray full_solution_arr

    # Love number information
    cdef cnp.ndarray complex_love_arr

    # EOS solution arrays
    cdef public cnp.ndarray gravity_array
    cdef public cnp.ndarray pressure_array
    cdef public cnp.ndarray mass_array
    cdef public cnp.ndarray moi_array
    cdef public cnp.ndarray density_array
    cdef public cnp.ndarray shear_modulus_array
    cdef public cnp.ndarray bulk_modulus_array

    cdef void set_model_names(
        self,
        int* bc_models_ptr) noexcept nogil
    cdef change_radius_array(
        self,
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        cpp_bool array_changed = *) noexcept
