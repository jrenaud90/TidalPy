from TidalPy.utilities.classes.base_x cimport TidalPyBaseExtensionClass


cdef class RheologyModelBase(TidalPyBaseExtensionClass):

    # Class attributes
    cdef size_t expected_num_args

    # Class methods
    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil

    cdef void _vectorize_frequency(
            self,
            double* frequency_ptr,
            double modulus,
            double viscosity,
            double complex* output_ptr,
            Py_ssize_t n,
            ) noexcept nogil

    cdef void _vectorize_modulus_viscosity(
            self,
            double frequency,
            double* modulus_ptr,
            double* viscosity_ptr,
            double complex* output_ptr,
            Py_ssize_t n,
            ) noexcept nogil
