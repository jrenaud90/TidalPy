import numpy as np
cimport numpy as np
np.import_array()


from libcpp cimport bool as cpp_bool

from CyRK.utils.vector cimport vector

cdef extern from "solutions_.cpp" nogil:
    
    const int MAX_NUM_Y
    const int MAX_NUM_Y_REAL

    cdef cppclass RadialSolutionStorageCC:
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

        vector[double] complex_love_vec
        double* complex_love_ptr
        vector[double] full_solution_vec
        double* full_solution_ptr

        void find_love(double surface_gravity) noexcept nogil
        set_message(const char* new_message) noexcept nogil


cdef class RadialSolverSolution:

    # Size and state information
    cdef size_t num_slices
    cdef char num_ytypes
    cdef cpp_bool ytype_names_set
    cdef char* ytypes[5]

    # Main storage container
    cdef RadialSolutionStorageCC* solution_storage_ptr

    # Result pointers and data
    cdef np.ndarray full_solution_arr

    # Love number information
    cdef np.ndarray complex_love_arr

    cdef void set_model_names(self, int* bc_models_ptr) noexcept nogil


cdef size_t cf_find_num_shooting_solutions(
    int layer_type,
    bint is_static,
    bint is_incompressible
    ) noexcept nogil
