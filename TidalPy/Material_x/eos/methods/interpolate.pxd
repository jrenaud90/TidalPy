from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Material_x.eos.ode cimport c_EOS_ODEInput


cdef extern from "interpolate_.hpp" nogil:

    cdef struct c_InterpolateEOSInput:
        size_t  num_slices
        double* radius_array_ptr
        double* density_array_ptr
        cpp_complex[double]* bulk_modulus_array_ptr
        cpp_complex[double]* shear_modulus_array_ptr


    cdef void c_preeval_interpolate(
            # Values that will be updated by the function
            char* preeval_output,
            # Input that is used by the pre-eval
            double radius,
            double* radial_solutions,
            char* preeval_input
            ) noexcept nogil
