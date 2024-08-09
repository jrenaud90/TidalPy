

import numpy as np
cimport numpy as np
np.import_array()

from CyRK.utils.memory cimport shared_ptr, make_shared

# Base Class
cdef extern from "base_.cpp" nogil:
    cppclass RheologyModelBaseCC:
        char rheology_model_int
        char num_parameters
        char* rheology_name_ptr
    
        RheologyModelBaseCC()
        RheologyModelBaseCC(
            const double* parameters,
            int num_parameters, 
            char rheology_model_int,
            char* rheology_name)
        void update_parameters(
            const double* parameters)
        void call(
            double* output_ptr,
            double frequency,
            double modulus,
            double viscosity)
        void vectorize_frequency(
            double* output_ptr,
            size_t output_size,
            double* frequency_ptr,
            double modulus, 
            double viscosity)
        void vectorize_modulus_viscosity(
            double* output_ptr,
            size_t output_size,
            double frequency,
            double* modulus_ptr, 
            double* viscosity_ptr)
        void vectorize_all(
            double* output_ptr,
            size_t output_size,
            double* frequency_ptr,
            double* modulus_ptr, 
            double* viscosity_ptr)


# Rheology Model (subclass)
cdef extern from "models_.hpp" nogil:
    const char ELASTIC_RHEOLOGY_INT
    const char ELASTIC_NUM_PARAMETERS
    cppclass Elastic(RheologyModelBaseCC):
        pass
    
    const char NEWTON_RHEOLOGY_INT
    const char NEWTON_NUM_PARAMETERS
    cppclass Newton(RheologyModelBaseCC):
        pass

    const char MAXWELL_RHEOLOGY_INT
    const char MAXWELL_NUM_PARAMETERS
    cppclass Maxwell(RheologyModelBaseCC):
        pass

    const char VOIGT_RHEOLOGY_INT
    const char VOIGT_NUM_PARAMETERS
    cppclass Voigt(RheologyModelBaseCC):
        pass

    const char BURGERS_RHEOLOGY_INT
    const char BURGERS_NUM_PARAMETERS
    cppclass Burgers(RheologyModelBaseCC):
        pass
    
    const char ANDRADE_RHEOLOGY_INT
    const char ANDRADE_NUM_PARAMETERS
    cppclass Andrade(RheologyModelBaseCC):
        pass
    
    const char SUNDBERG_RHEOLOGY_INT
    const char SUNDBERG_NUM_PARAMETERS
    cppclass SundbergCooper(RheologyModelBaseCC):
        pass


cdef class RheologyModel:

    cdef char rheology_int
    cdef char required_num_params
    cdef str rheology_name
    cdef RheologyModelBaseCC* rheology_class_ptr
    

cdef shared_ptr[RheologyModelBaseCC] cf_find_and_build_rheology(
    const char rheology_int,
    double* rheology_parameters) noexcept nogil

cdef char cf_get_rheology_int(
    const char* rheology_name) noexcept nogil