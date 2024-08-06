from libcpp cimport bool as cpp_bool


from CyRK.utils.memory cimport shared_ptr

from TidalPy.rheology.models cimport RheologyModelBaseCC

cdef struct interpolate_input:
    size_t num_slices
    double* radius_array_ptr
    double* density_array_ptr
    double* bulk_modulus_array_ptr
    double* bulk_viscosity_array_ptr
    double* shear_modulus_array_ptr
    double* shear_viscosity_array_ptr
    shared_ptr[RheologyModelBaseCC] bulk_rheology
    shared_ptr[RheologyModelBaseCC] shear_rheology

cdef void preeval_interpolate(
        # Values that will be updated by the function
        void* preeval_output,
        # Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        void* preeval_input
        ) noexcept nogil