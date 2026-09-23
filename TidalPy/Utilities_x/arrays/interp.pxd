# distutils: language = c++

from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector


cdef extern from "interp_.hpp" namespace "tidalpy" nogil:
    double c_interp(
        double desired_x,
        const double* x_domain,
        const double* dependent_values,
        size_t len_x,
        size_t guess)

    cpp_complex[double] c_interp_complex(
        double desired_x,
        const double* x_domain,
        const cpp_complex[double]* dependent_values,
        size_t len_x,
        size_t guess)


cdef extern from "layer_partition_.hpp" namespace "tidalpy" nogil:
    void c_partition_radius_by_layer(
        const double* radius_ptr,
        size_t num_slices,
        const double* upper_radius_bylayer_ptr,
        size_t num_layers,
        vector[size_t]& first_slice_out,
        vector[size_t]& num_slices_out)
